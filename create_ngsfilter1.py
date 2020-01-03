import pandas as pd
import os
import sys
from itertools import repeat
import click

# Create NGS filter
# Glossary:
# PP = primer plate, combination of tags and primers. Tags are unique for each well. Eight plates per library.
# AP = aliquot plate, each well holds its own sample.
# IMPORTANT: AP files must be named as follow: *****_LIBRARYNAME_AP.xls(x).
# Libraries, APs and PPs must have the same names in all input files.
@click.command()
@click.option("--project", help="Path to project folder")
@click.option("--platenames", help="Filename for excel table of aliquot plate and primer plates mapping")
@click.option("--tagnames", default="Tagscombo_UA.csv", help="Tag combination, filename")
@click.option("--primernames", default="primer_list.csv", help="Primers for each loci, filename")

def main(project, platenames, tagnames, primernames):
    if not os.path.exists(project):
        print("Incorrect project folder. Please, check your path!")
        exit()
    os.chdir(project)
    output_path = "./1_ngsfilters/"
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    if not os.path.exists("./0_prep_ngsfilters/" + platenames):
        print("Aliquot plate file is not found (--platenames argument is incorrect)."
              "File must be located in project/0_prep_ngsfilters/. You don't need provide a path, only a filename.")
        exit()
    AP_map = pd.read_excel("./0_prep_ngsfilters/" + platenames)

    if not os.path.isfile("./0_prep_ngsfilters/" + tagnames):
        print("Tag combination file is not found (--tagnames argument is incorrect).\n"
              "File must be located in project/0_prep_ngsfilters/. You don't need provide a path, only a filename.")
        exit()
    tagcombo = pd.read_csv("./0_prep_ngsfilters/" + tagnames)

    if not os.path.isfile("./0_prep_ngsfilters/" + primernames):
        print("Primers file is not found (--primers argument is incorrect).\n"
              "File must be located in project/0_prep_ngsfilters/. You don't need provide a path, only a filename.")
        exit()
    primers = pd.read_csv("./0_prep_ngsfilters/" + primernames)

    AP_path = "./0_prep_ngsfilters/aliquot_plates/"

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
                    self.APfile = pd.read_excel(AP_path + self.AP)
                    self.APfile["tags"] = tagcombo[self.PP]
                    self.APfile["SPositionId"] = [x for x in range(1,97)] # can we do in for non-96 wells??
                    self.intermediate = self.APfile.apply(lambda row: pd.DataFrame(list(zip(primers.iloc[:, 0], repeat(row["SPositionBC"]),
                                                           repeat(row["SPositionId"]), repeat(self.PP), repeat(row["tags"]),
                                                           primers.iloc[:, 1], primers.iloc[:, 2], repeat("F"),
                                                           repeat("@")))), axis=1)
                    self.ngsfilters.append(pd.concat(self.intermediate.tolist()))
            self.ngsfilter = pd.concat(self.ngsfilters)
            self.ngsfilter[2] =  self.ngsfilter[2].apply(lambda x: "%04d" % x) # make position standard - do we really need this??
            self.ngsfilter[1] = self.ngsfilter[[1, 2, 3]].apply(lambda row: '_'.join(row.values.astype(str)), axis=1) # prepare sample names in format ID_position_PP
            self.ngsfilter = self.ngsfilter.sort_values([3, 2, 0]) # sort by PP, position and loci name
            self.ngsfilter = self.ngsfilter.drop(columns=[2, 3]) # remove extra columns
        def write_ngsfilter(self):
            self.ngsfilter.to_csv(output_path + self.name + ".ngsfilter", sep="\t", header=False, index=False)

    libraries = AP_map["Library_BC"].unique() # list of libraries

    for lib in libraries:
        run = Library(lib)
        print(f"{lib} initialized.")
        run.create_ngsfilter()
        print(f"Creating ngsfilter for {lib}.")
        run.write_ngsfilter()

if __name__ == '__main__':
    main()