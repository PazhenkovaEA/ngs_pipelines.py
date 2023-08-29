import pandas as pd
import os
import sys
from itertools import repeat
import click

# Create NGS filter
# Glossary:
# PP = primer plate, combination of tags and primers. Tags are unique for each well. Eight plates per library.
# AP = aliquot plate, each well holds its own sample.
# IMPORTANT: AP files must be named as follow: AP_LIBRARYNAME_APName.xls(x).
# Libraries, APs and PPs must have the same names in all input files.
@click.command()
@click.option("--project", help="Path to project folder")
@click.option("--plates", help="Path to the excel table of aliquot plate and primer plates mapping")
@click.option("--tags", help="Tag combination, path to the .csv file")
@click.option("--primers", help="Primers for each loci, path to the .csv file")
@click.option("--aliquotplates", default="AP", help="Folder name (not a path) with Aliquot plates, should be placed to the Project folder")

def main(project, plates, tags, primers, aliquotplates):
    if not os.path.exists(project):
        print("Incorrect project folder. Please, check your path!")
        exit()
    #os.chdir(project)
    output_path = f"{project}/ngsfilters/"
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    result=f"{project}/results"
    if not os.path.exists(result):
        os.mkdir(result)

    if not os.path.isfile(plates):
        print("Aliquot plate file is not found (--plates argument is incorrect). Please, check your path")
        exit()
    AP_map = pd.read_excel(plates)

    if not os.path.isfile(tags):
        print("Tag combination file is not found (--tags argument is incorrect). Please, check your path")
        exit()
    tagcombo = pd.read_csv(tags)

    if not os.path.isfile(primers):
        print("Primers file is not found (--primers argument is incorrect). Please, check your path")
        exit()
    primers = pd.read_csv(primers)
    primers = primers[["locus",	"primerF", "primerR"]].drop_duplicates()

    AP_path = f"{project}/{aliquotplates}"
    if not os.path.exists(AP_path):
        print("Folder with aliquot plates is not found (--aliquotplates argument is incorrect).")
        exit()

    aliquot_plates = sorted(os.listdir(AP_path))
    if len(aliquot_plates) == 0:
        print("Aliquot plates are not found. Please, check your paths!")
        sys.exit()

    class Library:
        def __init__(self, libname=None):
            self.name = libname
            self.aliquot_plates = [plate for plate in aliquot_plates if self.name in plate] # filter aliquot plates for current library
            # this dictionary is just about mapping of PP to corresponding aliquot plate file. Here be careful with AP filenames and with column names in resulted file for mapping PP to AP.
            self.AP_to_PP = {x:sorted(AP_map[(AP_map["Library_BC"] == x.split(".")[0].split("_")[-2])
                                             & (AP_map["Aliquot Plate"] == x.split(".")[0].split("_")[-1])]["Primer Plate"].values.tolist()) for x in self.aliquot_plates}
            self.ngsfilters = []
            self.ngsfilter = None

        def create_ngsfilter(self):
            for self.AP in self.AP_to_PP.keys():
                for self.PP in self.AP_to_PP[self.AP]:
                    self.APfile = pd.read_excel(AP_path + "/" + self.AP)
                    self.APfile["tags"] = tagcombo[self.PP]
                    #self.APfile["SPositionId"] = [x for x in range(1,97)] # can we do it for non-96 wells??
                    self.intermediate = self.APfile.apply(lambda row: pd.DataFrame(list(zip(primers.iloc[:, 0], repeat(row["SPositionBC"]),
                                                           repeat(row["SPositionId"]), repeat(self.PP), repeat(row["tags"]),
                                                           primers.iloc[:, 1], primers.iloc[:, 2], repeat("F"),
                                                           repeat("@")))), axis=1)
                    self.ngsfilters.append(pd.concat(self.intermediate.tolist()))
            self.ngsfilter = pd.concat(self.ngsfilters)

            self.position = self.ngsfilter[[1,3,0,2,4]].copy() #Create a file with sample position.
            self.position.columns = ['Sample_Name', 'Plate', 'Marker','Position', 'TagCombo']
            self.position["Plate"] = self.position["Plate"].apply(lambda row: row.replace("PP", ""))
            self.position["Run_Name"] = self.name
            self.position['length'] = ""
            self.position['Read_Count'] = ""
            self.listcol = ['Sample_Name', 'Plate', 'Read_Count', 'Marker', 'Run_Name', 'length', 'Position', 'TagCombo']
            self.position = self.position.reindex(self.listcol, axis=1)
            self.ngsfilter[1] = self.ngsfilter[[1, 2, 3]].apply(lambda row: '__'.join(row.values.astype(str)), axis=1) # prepare sample names in format ID_position_PP
            self.ngsfilter = self.ngsfilter.sort_values([3, 2, 0]) # sort by PP, position and loci name
            self.ngsfilter = self.ngsfilter.drop(columns=[2, 3]) # remove extra columns
        def write_ngsfilter(self):
            self.ngsfilter.to_csv(output_path + self.name + ".ngsfilter", sep="\t", header=False, index=False)
        def write_positions(self):
            self.position.to_csv(f"{result}/{lib}_sample_positions.txt", mode="a", sep="\t", index=False)



    libraries = AP_map["Library_BC"].unique() # list of libraries

    for lib in libraries:
        run = Library(lib)
        print(f"{lib} initialized.")
        run.create_ngsfilter()
        print(f"Creating ngsfilter for {lib}.")
        run.write_ngsfilter()
        run.write_positions()

if __name__ == '__main__':
    main()



