#Graphs NGS
#Set by user
library_path = "/Users/elena/PycharmProjects/ngs_pipelines/LWA_workshop/DAB083"
proj = "DAB083"
plates = c(1:8) #primer plates numbers


library(ggplot2)
library(dplyr)

fpath = paste(library_path, "results/GenotypeExportTemp.txt", sep = "/")
data.types = c("character","character","character","character","character","character","character","character","numeric","character","character","character","character","numeric")

genotypeData=read.table(fpath, sep="\t", header=T, colClasses=data.types)
names(genotypeData) = c("Allele1","Allele2","Sample_Name","Marker","Plate","Allele","called","flag","Read_Count","stutter","Run_Name","Position","TagCombo","Length")
genotypeData$Sample_Name <- apply( genotypeData[ , c("Sample_Name", "Position")] , 1 , paste , collapse = "_" )
genotypeData$Sample_Name <- as.factor(genotypeData$Sample_Name)

genotypeData$IsFlagged = genotypeData$flag != ""
genotypeData$Titles = as.factor(paste("L:",genotypeData$Marker," *** G: ",genotypeData$Allele1," -- ",genotypeData$Allele2, sep=""))

#added
genotypeData$Length[genotypeData$Length == 0] = NA

genotypeData$Marker = as.factor(genotypeData$Marker)

genotypeData$flag[genotypeData$IsFlagged] = paste("*", genotypeData$flag[genotypeData$IsFlagged], sep="")

#mywait()

PlateID = 4 #number of primer plate to highlight
imagePath = paste(library_path,"results/genotypes.pdf", sep ="/")
pdf(imagePath, width = 50, height = 40)
for (level in unique(genotypeData$Sample_Name)) {
  # Subset the data for the current factor level
  subset_data <- subset(genotypeData, Sample_Name == level)
  current_plot <- ggplot(subset_data, aes(y=Read_Count, x=Length)) +
    scale_fill_manual(values=c(NA,"green","red","blue"))+
    #geom_hline(aes(yintercept=100, linetype="dashed"),color="blue", alpha=0.8)+
    geom_line(aes(group=as.factor(TagCombo), linetype = as.factor(TagCombo) ),color="orange") + #, alpha=0.5)+
    geom_line(data=filter(subset_data, Plate == PlateID), color="red") + #extra line for plate
    scale_size_manual(values=c(1.5,2.5))+
    geom_point(aes(shape=IsFlagged, color=IsFlagged, fill=IsFlagged, size=IsFlagged)) + 
    scale_shape_manual(values=21:25) +
    #scale_fill_brewer(type="qual")+
    scale_color_manual(values=c("red","black","blue","magenta2"))+
    #scale_color_manual(values=c("black", "magenta2", "blue", "red", "purple","green","cyan","orange")) +
    geom_text(aes(label=paste(Allele, flag, sep=""), color=called, y =Read_Count, size=IsFlagged), size=3,nudge_x=-0.07, check_overlap=T, hjust="right", vjust = "bottom", alpha=0.6 ) +
    #geom_text(aes(label=paste(Allele1,Allele2, sep =" *** ")),x=-Inf,y=Inf,hjust=0,vjust=1,color="red",size=7,alpha=0.15) +
    #geom_text(aes(label=Allele2),x=Inf,y=-Inf,hjust=1,vjust=0,color="red",size=6,alpha=0.15) +
    facet_wrap(~Titles, ncol=5, scales="free") + theme_bw() +#scales="free_y") +
    theme(legend.position="none")+xlab("")+ylab("")+
    theme(strip.background =element_rect(fill="white"))+
    theme(strip.text = element_text(colour = 'red', size=11)) + 
    ggtitle(paste("Sample name:", level))
  print(current_plot)
}
dev.off()


#Read count plot
library(data.table)
library(ggplot2)
library(gridExtra)
showSumsByLibrary <- function(mc, loc = "./", pos = c("by_plate", "by_sample")) {
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  # change function sum to mean, if you're interested in mean
  mcp <- mc[, .(mean.count = sum(Read_Count, na.rm = TRUE)), by = .(Sample_Name, Run_Name, Plate, Position)]
  
  out <- sapply(split(mcp, f = list(mcp$Run_Name)), FUN = function(x) {
    
    xy <- split(mcp, f = mcp$Plate)
    
    # In this chunk we check the existance of every plate. If one is missing, we
    # put an NA value in its stead.
    platenames <- sapply(xy, FUN = function(m) unique(m$Plate))
    names(xy) <- platenames
    allplates <- vector("list", length(pos))
    names(allplates) <- platenames
    
    for (i in pos) {
      index <- as.character(i)
      tmp <- xy[index]
      # If a plate is missing, add NA to the entry
      if (is.null(tmp)) {
        allplates[[index]]<- index
      } else {
        allplates[[index]] <- data.frame(tmp)
        colnames(allplates[[index]]) <- names(tmp[[index]])
      }
      rm(tmp)  # this will make sure things are clean once loop is finished
    }
    
    out <- sapply(allplates, FUN = function(y) {
      if (length(y) == 1 & is.character(y)) {
        # If there is no data for a certain plate, add blank plate.
        blankpage <- ggplot() +
          geom_text(aes(x = 0.5, y = 0.5), label = "No data", size = 24) +
          theme_bw() +
          ggtitle(y) +
          geom_blank() +
          scale_x_continuous(breaks = seq(0, 1, length.out = 12), labels = 1:12) +
          scale_y_continuous(breaks = seq(0, 1, length.out = 8), labels = rev(LETTERS[1:8])) +
          theme(axis.title = element_blank(), panel.grid.minor = element_blank())
        
        return(blankpage)
      }
      # create a vector and add data for correct positions
      # this is necessary to preserve missing positions
      y$xpos <- as.numeric(y$Position)
      xout <- rep(NA, 96)
      lbs <- xout # vector used to fill in labels
      nms <- xout # vector used to fill in names
      xout[y$xpos] <- y[, "mean.count"]
      lbs[y$xpos] <- sprintf("%s\n%s", y[, "Sample_Name"], round(y[, "mean.count"], 0))
      lbs <- matrix(lbs, nrow = 8)
      
      nms[y$xpos] <- y[, "Sample_Name"]
      nms <- matrix(nms, nrow = 8)
      colnames(nms) <- 1:12
      rownames(nms) <- LETTERS[1:8]
      nms <- reshape2::melt(nms)
      nms$isnegcont <- NA
      sn <- y[grepl("(^AC\\..*$|^AK\\..*$|NeKo.*$)", "Sample_Name"), "Sample_Name"]
      nms[nms$value %in% sn, "isnegcont"] <- TRUE
      
      out <- matrix(xout, nrow = 8)
      colnames(out) <- 1:12
      rownames(out) <- LETTERS[1:8]
      out <- reshape2::melt(out)
      out$Var2 <- as.factor(out$Var2)
      out$Var1 <- factor(out$Var1, levels = rev(levels(out$Var1)), ordered = TRUE)
      
      rownames(lbs) <- LETTERS[1:8]
      colnames(lbs) <- 1:12
      lbs <- melt(lbs, value.name = "lbl")
      lbs$Var2 <- as.factor(lbs$Var2)
      lbs$Var1 <- factor(lbs$Var1, levels = rev(levels(lbs$Var1)), ordered = TRUE)
      
      out <- Reduce(function(x, y) merge(x, y, by = c("Var1", "Var2")), list(out, nms, lbs))
      
      nk <- sapply(na.omit(out[out$isnegcont == TRUE, c("Var1", "Var2")]), as.numeric, simplify = FALSE)
      nk <- as.data.frame(nk)
      
      ggplot(out, aes(x = Var2, y = Var1)) +
        theme_bw() +
        scale_x_discrete(position = "top") +
        theme(axis.ticks = element_blank(), legend.position = "none") +
        xlab("") + ylab("") +
        ggtitle(unique(y$Plate)) +
        geom_tile(aes(fill = value.x)) +
        geom_rect(data = nk, aes(xmin = Var2 - 0.5, xmax = Var2 + 0.5, ymin = Var1 - 0.5, ymax = Var1 + 0.5), 
                  fill = NA, color = "red", size = 1) +
        geom_text(aes(label = lbl), size = 3.5, color = "white")
      
    }, simplify = FALSE)
    
    # decide if plates will be ordered to check plate performance or sample performance
    if (is.numeric(pos)) {
      # stop if pos has insufficient number of positions specified
      stopifnot(length(pos) == length(out))
    }
    
    # plot to file, loc should have a trailing slash included
    if (!grepl(".*/$", loc)) {
      stop("Trailing slash in `pos` not found. Please check your path.")
    }
    
    pdf(file = sprintf("%splatecount_%s.pdf", loc, unique(x$Run_Name)), width = 20, height = 20)
    do.call(grid.arrange, c(grobs = out, nrow = 4, as.table = FALSE))
    dev.off()
  }, simplify = FALSE)
  return(out)
}


genotypes <- sprintf("./%s/results/%s_genotypes.txt", library_path, proj)  
xy1 <- fread(genotypes, stringsAsFactors = FALSE,
             colClasses = list(character = c(1, 4, 5, 7, 9, 11, 12),
                               numeric = c(2, 3, 6)))

showSumsByLibrary(xy1, loc = paste(library_path, 'results', sep = "/"), pos = plates)
