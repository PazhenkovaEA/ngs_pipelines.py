import pandas as pd
import numpy as np
import ast
import re
import os
import click


@click.command()
@click.option("--project", help="Path to project folder")
@click.option("--motifs", help="Path to file for repeats motifs for each loci", default=None)
@click.option("--snp", help="Path to file for SNP sequences", default=None)
@click.option("--clean", default=True, help="Remove non-called alleles")


def main(project, motifs, snp, clean):
    project = "/Users/elena/PycharmProjects/ngs_pipelines/LF"
    motifs = "/Users/elena/PycharmProjects/ngs_pipelines/motifs.csv"
    snp = "/Users/elena/PycharmProjects/ngs_pipelines/LF/SNP.csv"
    stutter = 0.18 # "Read proportion for stutter"
    disbalance = 0.33 # Disbalanced allele
    LowCount = 100  # Low reads number threshold
    AlleleWithNoStutterHeight = 100 # Allele without stutter read number threshold

   # Read table with motif data and add columns with necessary data for allele calling if absent.
    if motifs:
        if not os.path.isfile(motifs):
            print("Check the path to the microsattelite motifs file")
            exit()
        motifs = pd.read_csv(motifs, converters={'locus': lambda x: str(x)}) #user input
        if not "Stutter" in motifs.columns:
            motifs["Stutter"] = float(stutter)
        if not "Disbalance" in motifs.columns:
            motifs["Disbalance"] = float(disbalance)
        if not "LowCount" in motifs.columns:
            motifs["LowCount"] = int(LowCount)
        if not "AlleleWithNoStutterHeight" in motifs.columns:
            motifs["AlleleWithNoStutterHeight"] = int(AlleleWithNoStutterHeight)
        motifs.columns = ["Marker", "Repeat","Stutter","Disbalance","LowCount", "AlleleWithNoStutterHeight"] #mt table from fishbone package
        motifs["Repeat"] = motifs["Repeat"].apply(lambda row: row.lower())


    # TODO: make possible to use fasta input
    if snp:
        if not os.path.isfile(snp):
            print("Check the path to the snp file")
            exit()
        snp = pd.read_csv(snp, header=None) #input
        snp.rename(columns={0: "Marker", 1: "Seq"}, inplace=True)
        snp["SNP"] = snp.apply(lambda row: row["Marker"].split("-")[0], axis=1)
        snp["Seq"] = snp["Seq"].apply(lambda row: row.lower())
    countTS = 20  # reads that do not match exact sequence should have this number of reads,
    # otherwise they are discarded
    lowCount = 100  # if unrecognized as sex sequence and below this threshold, flag as "L"
    junkTS = 0.15  # when scanning for junk, if a sequence is less than this of max sequence, discard
    dbTS = 0.5  # if sequence is below dbTS relative to max sequence, flag as disbalanced
    TS = 0.05  # 5% of divergence between sequences is possible if read count < 100

    microsat = False
    variant = False

    # Each row in .tab file contains data for several samples.
    # This function parse it to human-readable format, when one row = one sample. Will be applied to each row in .tab file.
    def tab_to_genotypes(row):
        genotypes1 = pd.DataFrame.from_dict(row[3], orient='index')
        genotypes1["Sample_Name"] = genotypes1.index
        genotypes1["length"],genotypes1["Sequence"] = row[5],row[0]
        return genotypes1

    ### Fishbone functions ###

    def find_stutter(row):
        stutters = data[data["length"] == row["length"] - len(motif.unique()[0])]
        #if len(stutters) == 0:
            #return None #no stutters!
        found_stutters = None
        # Find all motif occurences in the sequence.
        coordinates =[]
        for mo in motif:
            coordinates += [m.start() for m in re.finditer(mo, row["Sequence"])]
        #Shorted candidate sequence by one motif
        for c in coordinates:
            sequence = row["Sequence"][0:c] + row["Sequence"][c+len(motif.unique()[0]):] #remove one motif
            if len(stutters[stutters["Sequence"] == sequence]) != 0:
                found_stutters = stutters[stutters["Sequence"] == sequence]
            else:
                continue
        return found_stutters


    def call_allele(row, index):
        rh = row["Read_Count"]/maxA
        # 1. if allele has number of reads < L, flag it as "L"
        # (performed towards the end)
        # 2. see if allele has a stutter
        stutter = find_stutter(row)
        # 2a. if yes, mark as called
        # pandas needs to make a copy of data to change it
        if stutter is not None:
            if (len(stutter) == 1) & (rh > S):
                alleles.loc[index, "called"] = True
                alleles.loc[stutter.index, 'stutter'] = True  # 2ab. mark stutter as such
            # 2aa. if A in disbalance (A < D), flag as "D"
                if rh < D:
                    alleles.loc[index, "flag"] = "".join([alleles.loc[index, "flag"], "D"])  # 9  - flag column
            else:
                if row["Read_Count"] >= N:
                    alleles.loc[index, "flag"] = "".join([alleles.loc[index, "flag"], "N"])  # 9  - flag column
        else:
        # 2b. if no, check AlleleWithNoStutterHeight
        # 2ba. if x > AlleleWithNoStutterHeight, add flag "N"
            if row["Read_Count"] >= N:
                alleles.loc[index, "flag"] = "".join([alleles.loc[index, "flag"], "N"]) # 9  - flag column
            else:
                pass

    # Needleman-Wunsh Algorithm-like method to score a difference between the sequences - for SNP-calling
    def alignment_score(seq1, seq2):
        n = len(seq1)
        m = len(seq2)
        Wmin = np.zeros((n + 1, m + 1), dtype=int)
        Wmin[0] = [i for i in range(m + 1)]
        Wmin[:, 0] = [i for i in range(n + 1)]
        for i in range(1, (n + 1)):
            for j in range(1, (m + 1)):
                a1 = Wmin[i, j - 1] + 1
                a2 = Wmin[i - 1, j] + 1
                a3 = Wmin[i - 1, j - 1] + int(seq1[i - 1] != seq2[j - 1])
                Wmin[i][j] = min([a1, a2, a3])
        return Wmin[n][m]

    def call_snp(row, index):
        test = row["Sequence"]
        reference["Score"] = reference.apply(lambda line: alignment_score(line["Seq"], test), axis=1)
        lvls = row["Read_Count"] / max(data["Read_Count"])
        if lvls <= junkTS:
            pass
        if min(reference["Score"]) == 0:
            alleles.loc[index, "called"] = True
        elif (min(reference["Score"]) <= round(len(reference["Seq"].unique()[0]) * TS)) & (row["Read_Count"] > countTS):
            alleles.loc[index, "called"] = True
        elif row["Read_Count"] > lowCount:
            alleles.loc[index, "called"] = True
        if lvls < dbTS:
            alleles.loc[index, "flag"] = "".join([alleles.loc[index, "flag"], "D"])

    tab_dirs = list(filter(lambda x:'tab' in x, os.listdir(project))) #extract folders with "tab" in name - the output of obitools3 part
    genotypes = [] # list with final genotypes
    # Run for each library
    for lib in tab_dirs:
        library = lib.split("_")[0]
        tabs = os.listdir(f"{project}/{lib}")
        # Check, if ngsfilter exists.
        ngsfilter_path = f"{project}/1_ngsfilters/{library}.ngsfilter"
        if not os.path.isfile(ngsfilter_path):
            print(f"{library}.ngsfilter file have not found. Please, place it to the {project}/1_ngsfilters directory and try again.")
            exit()
        for t in tabs:
            # 1. Reformat obitools output
            #locus = t.split(".")[0].split("_")[-1]  # locus names
            locus = t.split("_")[-1]  # locus names
            if not library in t:
                continue
            if motifs is not None:
                if any(motifs["Marker"].str.contains(locus)):  # skip nonexistent loci
                    microsat = True
            elif snp is not None:
                if any(snp["SNP"].str.contains(locus)):
                    variant = True
            else:
                continue
            print(f"Processing {locus} locus, {library}")

            tab = pd.read_csv(f"{project}/{lib}/{t}", sep = "\t", header=None) # Columns: 0 - sequence, 3 - dict with samples, 4 - total count, 5 - length
            tab[3] = tab[3].apply(lambda row: ast.literal_eval(row)) # Transform obitools output in format {samplename:read_count} to dictionary
            intermediate = tab.apply(lambda row: tab_to_genotypes(row), axis=1) # Creates dataframe-like object, which contains dataframes with genotypes
            gen = pd.concat(intermediate.tolist()) # Combine it to one dataframe.
            #Add TagCombo
            ngsfilter = pd.read_csv(f"{project}/1_ngsfilters/{library}.ngsfilter", sep="\t", header=None) #specify path as variable
            ngsfilter = ngsfilter[ngsfilter[0].apply(lambda row: row.split("_")[-1]) == locus] #select only == locus
            ngsfilter[["Sample","TagCombo"]] = ngsfilter[[1, 2]]
            ngsfilter = ngsfilter[["Sample","TagCombo"]]
            gen = gen.merge(ngsfilter, left_on="Sample_Name", right_on="Sample")
            gen = gen.drop(["Sample"], axis=1)
            gen.rename(columns={0:"Read_Count", 2:"TagCombo"}, inplace=True)
            gen["Plate"] = gen["Sample_Name"].apply(lambda row: row.split("_")[-1])
            gen["Position"] = gen["Sample_Name"].apply(lambda row: row.split("_")[-2])
            gen["Sample_Name"] = gen["Sample_Name"].apply(lambda row: str(row.split("_")[0]))  # PLEASE, don't use underscore _ in sample names
            gen["Marker"], gen["Run_Name"] = locus, library
            gen.sort_values(by="Read_Count", inplace=True, ascending=False)
            #print(gen["Read_Count"].sum())
            #gen = gen[gen["Read_Count"] >= 9]
            split = list(gen.groupby(["Sample_Name", 'Plate','Position']))  # full version is overkill["Sample_name", 'Marker', 'Plate_name', 'Library', 'Position']
                # Allele calling for each sample
            for data in split:
                data = data[1]
                data.index = range(0, len(data)) # Create numeric index
                alleles = data.copy() # Create output dataframe
                alleles["called"], alleles["flag"], alleles["stutter"] = False, "", False # Create columns with genotype info
                if microsat:
                    motif = motifs[motifs["Marker"] == locus]["Repeat"] # Motif, corresponding to selected locus MAKE IT "MOTIFS"!
                    S = float(motifs[motifs["Marker"] == locus]["Stutter"].unique()[0])  # stutter
                    D = float(motifs[motifs["Marker"] == locus]["Disbalance"].unique()[0])  # disbalance
                    N = int(motifs[motifs["Marker"] == locus]["AlleleWithNoStutterHeight"].unique()[0])
                    L = int(motifs[motifs["Marker"] == locus]["LowCount"].unique()[0])  # low amplification threshold
                    maxA = max(data["Read_Count"])
                    for index,row in alleles.iterrows():
                        call_allele(row, index)
                    # 3. if number of unflagged A > 2, add flag "M" to all alleles
                    if not alleles["called"].any():
                        continue # if there are no alleles, go forward
                    if len(alleles[alleles["called"] & (alleles["flag"] != "D")]) > 2:
                        alleles["flag"] = alleles.apply(
                            lambda x: "".join([x["flag"], "M"]) if (x["called"] == True) & (x["flag"] == "") else x["flag"],
                            axis=1)
                    # 1. if allele has number of reads < L, flag it as "L"
                    alleles["flag"] = alleles.apply(lambda x: "".join([x["flag"], "L"]) if (x["Read_Count"]< L) & (x["called"] == True) else x["flag"], axis=1)
                if variant:
                    reference = snp.loc[snp["SNP"] == locus].copy()
                    for index, row in alleles.iterrows():
                        call_snp(row, index)
                genotypes.append(alleles)
            microsat = False
            #variant = False

        all_geno = pd.concat(genotypes)
        all_geno.rename(columns={'Read_count':"Read_Count", 'Sample_name': 'Sample_Name','Length':"length", 'Plate_name':"Plate",'Library': "Run_Name"}, inplace=True)
        listcol = ["Sample_Name", "Plate",	"Read_Count",	"Marker",	"Run_Name",	"length",	"Position",	"called",	"flag",	"stutter",	"Sequence", "TagCombo"]
        # If clean == TRUE, return only sequences which were tagged as allele or stutter
        if clean:
            all_geno = all_geno[all_geno["stutter"] | all_geno["called"]]
        all_geno = all_geno.reindex(listcol, axis=1)
        all_geno["Plate"] = all_geno["Plate"].apply(lambda row: row.replace("PP", ""))
        all_geno["Position"] = all_geno["Position"].apply(lambda row: int(row))
        all_geno["Position"] = all_geno["Position"].apply(lambda row: str(row))
        all_geno["called"] = all_geno["called"].apply(lambda row: str(row).replace("alse", "ALSE"))
        all_geno["called"] = all_geno["called"].apply(lambda row: str(row).replace("rue", "RUE"))
        all_geno["stutter"] = all_geno["stutter"].apply(lambda row: str(row).replace("alse", "ALSE"))
        all_geno["stutter"] = all_geno["stutter"].apply(lambda row: str(row).replace("rue", "RUE"))
        all_geno.to_csv(f"{project}/{library}_genotypes.txt", sep="\t", index=False)

        frequency = all_geno.groupby(["Sequence", "Marker"])["Read_Count"].sum()
        frequency = frequency.to_frame().reset_index().sort_values(["Marker", "Read_Count"], ascending=[True, False]).rename(columns={"Read_Count":"N"}).reindex(["Marker", "N", "Sequence"])
        frequency.to_csv(f"{project}/{library}_frequency_of_sequences_by_marker.txt", sep="\t", index=False)

if __name__ == '__main__':
    main()