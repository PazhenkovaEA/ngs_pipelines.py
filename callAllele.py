import pandas as pd
import numpy as np
import ast
import re
import os
import click
import subprocess
from collections import Counter

@click.command()
@click.option("--project", help="Path to project folder")
@click.option("--primers", help="Path to file with repeats motifs for each loci", default=None)
@click.option("--clean", default=True, help="Remove non-called alleles")
@click.option("--progressive_threshold", default=False, help="Threshold will be set for each loci equal to the maximum negative control read count")
@click.option("--make_graphs", default=True, help="Create a pdf file with graphs to manually check consensus")
@click.option("--reference_alleles", default="", help="Path to reference allele names")
#set_value
def main(project, primers, clean, progressive_threshold,make_graphs, reference_alleles):
    #project = "/Users/elena/PycharmProjects/ngs_pipelines/DAB065" #user input
    #primers = "/Users/elena/PycharmProjects/ngs_pipelines/UA_primers.csv" #user input
    #progressive_threshold = False
    #clean = True
    stutter = 0.18 # "Read proportion for stutter"
    disbalance = 0.33 # Disbalanced allele
    LowCount = 100  # Low reads number threshold
    AlleleWithNoStutterHeight = 100 # Allele without stutter read number threshold

    result = f"{project}/results"
    if not os.path.exists(result):
        os.mkdir(result)


    # Read table with motif data and add columns with necessary data for allele calling if absent.
    if not os.path.isfile(primers):
        print("Check the path to the microsattelite motifs file")
        exit()

    if not os.path.isfile(reference_alleles):
        ref_allele = False
        print("WARNING: reference allele table is not found. Alleles will be named within this library. Reference table will be created.")
    else:
        ref_allele = True

    # Threshold for consensus
    AlleleAcceptanceThreshold = 2  # default by Tomaž
    AlleleAcceptanceThreshold_hetero = 2


    all_loci = pd.read_csv(primers, converters={'locus': lambda x: str(x)})
    motifs = all_loci.loc[all_loci["type"] == "microsat"].copy()
    if len(motifs) != 0:
        if not "Stutter" in motifs.columns:
            motifs["Stutter"] = float(stutter)
        if not "Disbalance" in motifs.columns:
            motifs["Disbalance"] = float(disbalance)
        if not "LowCount" in motifs.columns:
            motifs["LowCount"] = int(LowCount)
        if not "AlleleWithNoStutterHeight" in motifs.columns:
            motifs["AlleleWithNoStutterHeight"] = int(AlleleWithNoStutterHeight)
        #motifs.columns = ["locus", "Repeat","Stutter","Disbalance","LowCount", "AlleleWithNoStutterHeight"] #mt table from fishbone package
        motifs["motif"] = motifs["motif"].apply(lambda row: row.lower())
    else:
        motifs = None

    snp = all_loci.loc[all_loci["type"] == "snp"].copy()
    if len(snp) != 0:
        snp["sequence"] = snp["sequence"].apply(lambda row: row.lower())
        countTS = 20  # reads that do not match exact sequence should have this number of reads,
        # otherwise they are discarded
        lowCount = 100  # if unrecognized as sex sequence and below this threshold, flag as "L"
        junkTS = 0.15  # when scanning for junk, if a sequence is less than this of max sequence, discard
        dbTS = 0.5  # if sequence is below dbTS relative to max sequence, flag as disbalanced
        TS = 0.05  # 5% of divergence between sequences is possible if read count < 100
        T = 1
    else:
        snp = None

    microsat = False
    variant = False

    # Each row in .tab file contains data for several samples.
    # This function parse it to human-readable format, when one row = one sample. Will be applied to each row in .tab file.
    def tab_to_genotypes(row):
        genotypes1 = pd.DataFrame.from_dict(row[3], orient='index')
        genotypes1["Sample_Name"] = genotypes1.index
        genotypes1["length"], genotypes1["Sequence"] = row[5], row[0]
        return genotypes1

    ### Fishbone functions ###

    def find_stutter(row):
        stutters = data[data["length"] == row["length"] - len(motif.unique()[0])]
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
        # 1. see if allele has a stutter
        stutter = find_stutter(row)
        # 2. if A has one stutter and its read count proportion higher than 0.18 from max - its called
        if stutter is not None:
            if (len(stutter) == 1) & (rh > S):
                alleles.loc[index, "called"] = True
                alleles.loc[stutter.index, 'stutter'] = True  # 2ab. mark stutter as such
            # 2aa. if A in disbalance (A < D), flag as "D"
                if rh < D:
                    alleles.loc[index, "flag"] = "".join([alleles.loc[index, "flag"], "D"])  # 9  - flag column
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
        reference["Score"] = reference.apply(lambda line: alignment_score(line["sequence"], test), axis=1)
        lvls = row["Read_Count"] / max(data["Read_Count"])
        if lvls <= junkTS * max(data["Read_Count"]):
            pass
        if min(reference["Score"]) == 0:
            alleles.loc[index, "called"] = True
        elif (min(reference["Score"]) <= round(len(reference["sequence"].unique()[0]) * TS)) & (row["Read_Count"] > countTS):
            alleles.loc[index, "called"] = True
        elif row["Read_Count"] > lowCount:
            alleles.loc[index, "called"] = True
        if lvls < dbTS:
            alleles.loc[index, "flag"] = "".join([alleles.loc[index, "flag"], "D"])
            alleles.loc[index, "called"] = True

    tab_dirs = list(filter(lambda x:'tab' in x, os.listdir(project))) #extract folders with "tab" in name - the output of obitools3 part
    genotypes = [] # list with final genotypes
    # Run for each library
    for lib in tab_dirs:
        library = lib.split("_")[0]
        tabs = os.listdir(f"{project}/{lib}")
        # Check, if ngsfilter exists.
        ngsfilter_path = f"{project}/ngsfilters/{library}.ngsfilter"
        if not os.path.isfile(ngsfilter_path):
            print(f"{library}.ngsfilter file have not found. Please, place it to the {project}/1_ngsfilters directory and try again.")
            exit()
        for t in tabs:
            # 1. Reformat obitools output
            #locus = t.split(".")[0].split("_")[-1]  # locus names
            locus = t.split(".")[0].split("__")[-1]  # locus names
            if not library in t:
                continue
            if motifs is not None:
                if any(motifs["locus"].str.contains(locus)):  # skip nonexistent loci
                    microsat = True
                    motif = motifs[motifs["locus"] == locus]["motif"]  # Motif, corresponding to selected locus MAKE IT "MOTIFS"!
                    S = float(motifs[motifs["locus"] == locus]["Stutter"].unique()[0])  # stutter
                    D = float(motifs[motifs["locus"] == locus]["Disbalance"].unique()[0])  # disbalance
            if snp is not None:
                if any(snp["locus"].str.contains(locus)):
                    variant = True
            if not variant | microsat:
                continue
            print(f"Processing {locus} locus, {library}")

            tab = pd.read_csv(f"{project}/{lib}/{t}", sep = "\t", header=None) # Columns: 0 - sequence, 3 - dict with samples, 4 - total count, 5 - length
            tab[3] = tab[3].apply(lambda row: ast.literal_eval(row)) # Transform obitools output in format {samplename:read_count} to dictionary
            intermediate = tab.apply(lambda row: tab_to_genotypes(row), axis=1) # Creates dataframe-like object, which contains dataframes with genotypes
            gen = pd.concat(intermediate.tolist()) # Combine it to one dataframe.
            #Add TagCombo
            ngsfilter = pd.read_csv(f"{project}/ngsfilters/{library}.ngsfilter", sep="\t", header=None)
            ngsfilter = ngsfilter[ngsfilter[0] == locus] #select rows corresponding to locus
            ngsfilter[["Sample","TagCombo"]] = ngsfilter[[1, 2]]
            ngsfilter = ngsfilter[["Sample","TagCombo"]]
            gen = gen.merge(ngsfilter, left_on="Sample_Name", right_on="Sample")
            gen = gen.drop(["Sample"], axis=1)
            gen.rename(columns={0:"Read_Count", 2:"TagCombo"}, inplace=True)
            gen["Plate"] = gen["Sample_Name"].apply(lambda row: row.split("__")[-1])
            gen["Position"] = gen["Sample_Name"].apply(lambda row: row.split("__")[-2])
            gen["Sample_Name"] = gen["Sample_Name"].apply(lambda row: str(row.split("__")[0]))
            gen.sort_values(by="Read_Count", inplace=True, ascending=False)
            gen["Marker"] = locus
            gen["Library"] = library
            if microsat:
                if progressive_threshold:
                    try:
                        #L = 1.5 * int(max(gen[gen["Sample_Name"].str.contains("^B[0-9]{2}")]["Read_Count"])) #maximum number of reads in Blanks named for example B01, B02 etc
                        L = 1.5 * int(max(gen[gen["Sample_Name"].str.contains("Blank")]["Read_Count"]))
                        N = L
                    except:
                        L = int(motifs[motifs["locus"] == locus]["LowCount"].unique()[0])  # low amplification threshold
                        N = int(motifs[motifs["locus"] == locus]["AlleleWithNoStutterHeight"].unique()[0])
                else:
                    L = int(motifs[motifs["locus"] == locus]["LowCount"].unique()[0])  # low amplification threshold
                    N = int(motifs[motifs["locus"] == locus]["AlleleWithNoStutterHeight"].unique()[0])
            #print(gen["Read_Count"].sum())
            gen = gen[gen["Read_Count"] >= 5] #filter out samples with number of reads < 1 to save time
            split = list(gen.groupby(["Sample_Name", 'Plate','Position']))  # full version is overkill["Sample_name", 'Marker', 'Plate_name', 'Library', 'Position']
                # Allele calling for each sample
            for data in split:
                data = data[1]
                data.index = range(0, len(data)) # Create numeric index
                alleles = data.copy() # Create output dataframe
                alleles["called"], alleles["flag"], alleles["stutter"] = False, "", False # Create columns with genotype info
                if microsat:
                    maxA = max(data["Read_Count"])
                    for index,row in alleles.iterrows():
                        call_allele(row, index)
                    if not alleles["called"].any():
                        continue # if there are no alleles, go forward
                    # 3. if allele has number of reads < L, flag it as "L"
                    alleles["flag"] = alleles.apply(lambda x: "".join([x["flag"], "L"]) if (x["Read_Count"] < L) & (x["called"]== True) else x[
                            "flag"], axis=1)
                    # 4. if number of unflagged A > 2, add flag "M" to all alleles
                    if len(alleles[(alleles["called"]==True) & (alleles["flag"] == "")]) > 2:
                        alleles["flag"] = alleles.apply(
                            lambda x: "".join([x["flag"], "M"]) if (x["called"]== True) & (x["flag"] == "") else x["flag"],
                            axis=1)
                    # 5. If there are only 2 alleles and and one of the alleles has L flag
                    #then if it's above D thresold, L flag will be removed.
                    if len(alleles[alleles["called"]== True]) == 2:
                        if (len(alleles[alleles["flag"] == "L"]) == 1) & (len(alleles[(alleles["flag"] == "") & (alleles["called"])]) == 1):
                            alleles["flag"] = ""
                    # 6. Alleles with no stutter can be called it they are higher than N threshold
                    alleles["flag"] = alleles.apply(
                        lambda x: "".join([x["flag"], "N"]) if (x["called"] != True) & (x["stutter"]!= True) & (
                                    x["Read_Count"] > N) else x["flag"], axis=1)
                    alleles[alleles["flag"] == "N"]["called"] = True
                if variant:
                    reference = snp.loc[snp["locus"] == locus].copy()
                    for index, row in alleles.iterrows():
                        call_snp(row, index)
                    alleles = alleles[alleles["Read_Count"] > T]
                if clean:
                    alleles = alleles[alleles["stutter"] | alleles["called"]]
                genotypes.append(alleles)
            microsat = False
            variant = False

        all_geno = pd.concat(genotypes)
        all_geno.rename(columns={'Read_count':"Read_Count", 'Sample_name': 'Sample_Name','Length':"length", 'Plate_name':"Plate",'Library': "Run_Name"}, inplace=True)
        listcol = ["Sample_Name", "Plate",	"Read_Count",	"Marker",	"Run_Name",	"length",	"Position",	"called",	"flag",	"stutter",	"Sequence", "TagCombo"]
        # If clean == TRUE, return only sequences which were tagged as allele or stutter


        if clean:
            all_geno = all_geno[all_geno["stutter"] | all_geno["called"]]
        all_geno = all_geno[listcol]
        all_geno["Plate"] = all_geno["Plate"].apply(lambda row: row.replace("PP", ""))
        all_geno["Position"] = all_geno["Position"].apply(lambda row: int(row))
        #all_geno["Position"] = all_geno["Position"].apply(lambda row: str(row))
        all_geno["called"] = all_geno["called"].apply(lambda row: str(row).replace("alse", "ALSE"))
        all_geno["called"] = all_geno["called"].apply(lambda row: str(row).replace("rue", "RUE"))
        all_geno["stutter"] = all_geno["stutter"].apply(lambda row: str(row).replace("alse", "ALSE"))
        all_geno["stutter"] = all_geno["stutter"].apply(lambda row: str(row).replace("rue", "RUE"))
        all_geno.to_csv(f"{result}/{library}_genotypes.txt", sep="\t", index=False)

        frequency = all_geno.groupby(["Sequence", "Marker"])["Read_Count"].sum()
        frequency = frequency.to_frame().reset_index().sort_values(["Marker", "Read_Count"], ascending=[True, False]).rename(columns={"Read_Count":"N"}).reindex(["Marker", "N", "Sequence"], axis=1)
        frequency.to_csv(f"{result}/{library}_frequency_of_sequences_by_marker.txt", sep="\t", index=False)

    # visualization of alleles for each sample and locus (call R script)
    if make_graphs:
        alleles_to_graphs = []
        for data in list(all_geno.groupby(["TagCombo", 'Marker'])):
            data = data[1]
            data = data.sort_values(["called", "Read_Count"], ascending=[False, False])
            data.index = range(0, len(data))  # Create numeric index
            # detect duplicated lengths
            duplicates = data.duplicated(subset='length', keep=False)
            suffixes = pd.Series('', index=data.index)
            for i, is_duplicate in enumerate(duplicates):
                if is_duplicate:
                    suffixes[i] = '_' + str(duplicates[:i].sum() + 1)
            data['Allele'] = data['length'].astype(str) + suffixes
            data.index = range(0, len(data))  # Create numeric index
            if len(data[data["called"] == 1]) > 1:
                Allele1 = data[data["called"] == 1].loc[0]["Allele"]
                Allele2 = data[data["called"] == 1].loc[1]["Allele"]
            else:
                Allele1 = data[data["called"] == 1].loc[0]["Allele"]
                Allele2 = ""
            data["Allele1"] = Allele1
            data["Allele2"] = Allele2
            data = data[["Allele1", "Allele2", "Sample_Name", "Marker", "Plate", "Allele", "called", "flag", "Read_Count",
                         "stutter", "Run_Name", "Position", "TagCombo", "length"]]
            alleles_to_graphs.append(data)
        alleles1 = pd.concat(alleles_to_graphs)
        alleles1.called = alleles1.called.replace({"TRUE": 1, "FALSE": 0})
        alleles1.stutter = alleles1.stutter.replace({"TRUE": 1, "FALSE": 0})
        alleles1.to_csv(f"{result}/GenotypeExportTemp.txt", sep="\t", index=False)
        subprocess.call(["Rscript", "PlotNGSGenotype.R", result])

    # Create a table allele names or read it from path. Then rename alleles.
    if ref_allele:
        ref_table = pd.read_table(reference_alleles)
        complete_reference = pd.merge(ref_table, frequency[["Sequence", "Marker", "N"]], on=["Sequence", "Marker"],
                                      how="outer")
        complete_reference.loc[complete_reference['N_x'].isnull(), 'N_x'] = complete_reference['N_y']
        complete_reference = complete_reference.rename(columns={'N_x': "N"}).drop(columns="N_y")
        complete_reference["Length"] = complete_reference.Sequence.apply(len)
        complete_reference = complete_reference.sort_values(["Marker", "AlleleName", "Length", "N"],
                                                            ascending=[True, True, True, False])
    else:
        complete_reference = frequency.copy()
        complete_reference["Length"] = complete_reference.Sequence.apply(len)
    complete_reference["Variant"] = complete_reference.groupby(["Marker", "Length"]).cumcount() + 1
    complete_reference["AlleleName"] = complete_reference.apply(lambda x: "_".join([str(x["Length"]), str(x["Variant"])]) if x["Variant"] > 1 else str(x["Length"]), axis=1)
    complete_reference.to_csv(f"{result}/reference_alleles_upd_{library}", sep="\t", index=False)
    a = pd.merge(complete_reference[["Sequence", "Marker", "AlleleName"]], all_geno, on=["Sequence", "Marker"])

    # Make a consensus
    cons = pd.DataFrame(columns=["Sample", "Mrkr", "Al1", "Al2", "Al3", "Al4", "NcnfA1", "NCnfA2", "ConfirmedAlleles",
                                 "UnconfirmedAlleles", "NAmp", "NAmpOK", "Success", "ADO", "ADORate", "QualityIndex"])

    # Estimate a number of replicates
    ngsfilter = pd.read_csv('/Users/elena/PycharmProjects/ngs_pipelines/DIVJA088/ngsfilters/DIVJA088.ngsfilter', sep="\t", header=None)
    ngsfilter[["Marker", "Sample_Name"]] = ngsfilter[[0, 1]]
    ngsfilter["Replicate"] = ngsfilter["Sample_Name"].apply(lambda row: row.split("__")[-1])
    ngsfilter["Sample_Name"] = ngsfilter["Sample_Name"].apply(lambda row: row.split("__")[0])
    ngsfilter = ngsfilter[["Sample_Name", "Replicate"]].drop_duplicates()
    replicates = ngsfilter.groupby(["Sample_Name"], as_index=False)["Replicate"].count()
    a = pd.merge(a, replicates, on = "Sample_Name") # calculate number of replicates for each sample


    for data in list(a.groupby(["Sample_Name", 'Marker'])):
        b = data[1]
        Sample = b["Sample_Name"].unique()[0]
        Marker = b["Marker"].unique()[0]
        NAmp = b["Replicate"].unique()[0]

        unconfirmed = []  # list for unconfirmed alleles
        confirmed = {}  # dict for flagged but confirmed alleles and its counts
        consensus = []  # list for consensus alleles

        alleles = Counter(b["AlleleName"])  # Create a dict with counts for each allele.
        # Sort - the most frequent alleles (by value first) and then (by key) allele names (shorter alleles and more frequent variants first).
        alleles = dict(sorted(alleles.items(), key=lambda item: (-item[1], item[0])))

        if len(alleles) > 1:
            Threshold = AlleleAcceptanceThreshold_hetero  # set up a threshold for heterozygotes if there are 2 or more possible alleles.
        else:
            Threshold = AlleleAcceptanceThreshold
        # Check each allele if it correspond to criteria
        for allele in list(alleles.keys()):
            data = b[b["AlleleName"] == allele]
            if len(data[pd.isna(data["flag"])]) == 0:  # if all alleles are flagged - all of them are unconfirmed
                unconfirmed.append(allele)
                continue
            else:
                if len(data[~pd.isna(data[
                                         "flag"])]) > 0:  # if some alleles are flagged but unflagged are present as well - add it to "Confirmed alleles".
                    confirmed.update({allele: len(data[~pd.isna(data["flag"])])})
                if len(data) > Threshold:
                    consensus.append(allele)
                else:
                    unconfirmed.append(allele)

        for i in range(4 - len(consensus)):  # what if length of consensus > 4?
            consensus += [""]
        # Check if there are confirmed alleles among the accepted consensus and how many replicates were flagged.
        NcnfA1, NcnfA2 = "", ""
        if len(confirmed) > 0:
            if consensus[0] in confirmed.keys():
                NcnfA1 = confirmed[consensus[0]]
            if consensus[1] in confirmed.keys():
                NcnfA2 = confirmed[consensus[1]]
        # Calculate number of successful amplifications and allelic dropouts (ADO). ADO - in how many replicated heterozygotes were not detected.
        # Amplification was not OK when its not heterozygote and contain second allele bellow threshold (how I understood from Tomaze's script).
        # I defined it differently - OK is when there are 2 alleles for heterozygotes and one for homozygotes.
        if len(confirmed) == 2:
            NAmpOK = NAmp - abs(alleles[list(confirmed.keys())[0]] - alleles[list(confirmed.keys())[1]])
            ADO = abs(alleles[list(confirmed.keys())[0]] - alleles[list(confirmed.keys())[1]])
            ADORate = ADO / NAmpOK
        elif len(confirmed) == 1:
            NAmpOK = alleles[list(confirmed.keys())[0]]
            ADO = 0
            ADORate = 0
        else:
            NAmpOK = 0
            ADO = 0
            ADORate = 0
        QualityIndex = NAmpOK / NAmp
        Success = NAmpOK / NAmp * 100

        cons = cons.append( {"Sample": Sample, "Mrkr": Marker, "Al1": consensus[0], "Al2": consensus[1], "Al3": consensus[2],
             "Al4": consensus[3], "NcnfA1": NcnfA1, "NCnfA2": NcnfA2,
             "ConfirmedAlleles": ';'.join(list(confirmed.keys())), "UnconfirmedAlleles": ';'.join(unconfirmed),
             "NAmp": NAmp, "NAmpOK": NAmpOK, "Success": Success,
             "ADO": ADO, "ADORate": ADORate, "QualityIndex": QualityIndex}, ignore_index=True)


if __name__ == '__main__':
    main()

