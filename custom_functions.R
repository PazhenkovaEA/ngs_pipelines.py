#' For each library will create 8 plates per page. Each plate has 8x12 positions with corresponding sample 
#' name and read count of all sequences.
#' Use parameter `loc` to output result to a specific folder. Do not forget a trailing slash. Defaults to `getwd()`.
#' 
#' Parameter pos will direct how plates are position in the output. `by_plate` it will preserve the order where
#' plates are 1-4 in the first column and 5-8 in the second column. `by_sample` will organize plates to have
#' plates with same replication together. `pos` can also be a numeric integer vector of the same length as there
#' are number of plates. This will designated in which order the plates are to be plotted. Keep in mind that 
#' the plotting order is column-wise.
#' 
showSumsByLibrary <- function(mc, loc = "./", pos = c("by_plate", "by_sample")) {
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  # change function sum to mean, if you're interested in mean
  mcp <- mc[, .(mean.count = sum(Read_Count, na.rm = TRUE)), by = .(Sample_Name, Run_Name, Plate, Position)]
  
  out <- sapply(split(mcp, f = list(mcp$Run_Name)), FUN = function(x) {
    xy <- split(x, f = x$Plate)
    
    # In this chunk we check the existance of every plate. If one is missing, we
    # put an NA value in its stead.
    platenames <- sapply(xy, FUN = function(m) unique(m$Plate))
    names(xy) <- platenames
    allplates <- vector("list", 8)
    names(allplates) <- platenames
    
    for (i in 1:8) {
      tmp <- xy[[as.character(i)]]
      # If a plate is missing, add NA to the entry
      if (is.null(tmp)) {
        allplates[[i]] <- as.character(i)
      } else {
        allplates[[i]] <- tmp
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
      
      xout[y$xpos] <- y[, mean.count]
      lbs[y$xpos] <- sprintf("%s\n%s", y[, Sample_Name], round(y[, mean.count], 0))
      lbs <- matrix(lbs, nrow = 8)
      
      nms[y$xpos] <- y[, Sample_Name]
      nms <- matrix(nms, nrow = 8)
      colnames(nms) <- 1:12
      rownames(nms) <- LETTERS[1:8]
      nms <- reshape2::melt(nms)
      nms$isnegcont <- NA
      sn <- y[grepl("(^AC\\..*$|^AK\\..*$|NeKo.*$)", Sample_Name), Sample_Name]
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
    if (pos == "by_sample") out <- out[c(1, 7, 3, 5, 2, 8, 4, 6)]
    if (is.numeric(pos)) {
      # stop if pos has insufficient number of positions specified
      stopifnot(length(pos) == length(out))
      out <- out[pos]
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

#' Clean ZF alleles.
#' @param fb Data as organized from the ngs pipeline.
#' @param ts Threshold value under which sequences are discarded relative to the sequence
#' with the highest number of reads in one sample * plate combination.

cleanZF <- function(fb, ts, db) {
  if (nrow(fb) == 0) return(NULL)
  # omit all sequences that do not reach ts% of the highest read.  
  lvls <- fb$Read_Count/max(fb$Read_Count)
  fb <- fb[lvls > ts, ]
  
  # if an allele is significantly below a db threshold, flag it as disbalanced
  lvls <- fb$Read_Count/max(fb$Read_Count)
  fb[lvls < db, ]$flag <- paste(fb[lvls < db, ]$flag, "D", sep = "")
  # fb[lvls < db, flag := paste(flag, "D", sep = "")]
  fb
}

#' Count number of repeats for samples * locus combinations.
countSampleLocusRepeats <- function(ngs) {
  it <- sapply(ngs, FUN = function(x) {
    ri <- fread(x, header = FALSE)
    lb <- gsub("\\.ngsfilter", "", basename(x))
    
    namerun <- data.table(Sample_Name = gsub("^(.*)(_\\d+_)(PP\\d+)$", "\\1", ri[, V2]),
                          Marker = gsub("^.*_([[:alnum:]]+)$", "\\1", ri[, V1]),
                          runname = gsub("^(.*)(_\\d+_)(PP\\d+)$", "\\3", ri[, V2]),
                          ident = ri[, V3])
    head(namerun)
    
    # namerun <- namerun[!duplicated(ident), ]
    
    namerun[, .(Run_Name = lb, .N), by = .(Sample_Name, Marker)]
  }, simplify = FALSE)
  it <- rbindlist(it)
}

#' Create a list of samples and their position in a library
#' @param x Full path to the .ngsfilter file.
sampleMetadata <- function(x) {
  lb <- fread(x, header = FALSE)
  ptn <- "^(.*)_(\\d+)_PP(\\d+)$"
  xy <- data.table(
    Sample_Name = gsub(ptn, "\\1", lb[, V2]),
    Plate = gsub(ptn, "\\3", lb[, V2]),
    Read_Count = NA,
    Marker = gsub("^.*_([[:alnum:]]{2})$", "\\1", lb[, V1]),
    Run_Name = gsub("^(.*[[:alnum:]])\\.ngsfilter$", "\\1", basename(x)),
    length = NA,
    Position = gsub(ptn, "\\2", lb[, V2]),
    TagCombo = lb[, V3]
  )
  xy
}

#' @param x A data.frame or data.table of called genotypes.
#' @param ngs A vector of paths to ngs filters.
#' @param sex Path to data about sex genotypes based on ZF.
#' 
findUnamplifiedSamples <- function(x, ngs, sex) {
  libs <- unique(x$Run_Name)
  zf <- fread(sex, header = TRUE,
              colClasses = list(character = c(1, 4, 5, 7, 9, 11, 12),
                                numeric = c(2, 3, 6)))
  
  sapply(libs, FUN = function(y, x, ngs, zf) {
    onelib <- x[x$Run_Name == y, ]
    zf <- zf[Run_Name == y, ]
    onengs <- ngs[ngs$lib == y, "path"]
    onengs <- fread(as.character(onengs), header = FALSE)
    
    onelib <- rbind(onelib, zf)
    
    # find all samples and construct a proper data.frame
    rgx <- "(^.*)_(\\d+)_PP(\\d+)$"
    ngx <- data.table(Sample_Name = gsub(rgx, "\\1", onengs$V2),
                      Plate = as.numeric(gsub(rgx, "\\3", onengs$V2)),
                      Position = gsub(rgx, "\\2", onengs$V2),
                      Marker = gsub("UA_MxRout1_(.*)$", "\\1", onengs$V1),
                      Run_Name = y,
                      TagCombo = onengs$V3,
                      stringsAsFactors = FALSE)
    
    # find all unique samples that have been amplified
    findups <- onelib[, .(Sample_Name, Plate, Position, Marker, Run_Name, TagCombo)]
    onelib <- onelib[!duplicated(findups), ]
    
    # append all relevant columns and extract only those that have no reads
    out <- merge(ngx, onelib, sort = FALSE, all = TRUE)
    out <- out[is.na(Read_Count), names(onelib), with = FALSE]
    out
    
  }, x = x, ngs = ngs, zf = zf, simplify = FALSE)
  
}
