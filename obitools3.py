import click
import pandas as pd
import os

@click.command()
@click.option("--project", help="Path to project folder")
@click.option("--library", help="Library name")
@click.option("--reads1", help="Forward .fastq reads")
@click.option("--reads2", help="Reverse .fastq reads")
@click.option("--primernames", default="primer_list.csv", help="Primers for each loci, filename")
@click.option("--obipath", default="~/obitools3", help="Path to obitools3 folder")


def main(project, library, obipath, reads1, reads2, primernames):
    os.chdir(project)
    if not os.path.exists(f"./{library}_tab"):
        os.mkdir(f"{library}_tab")
    if not os.path.isfile("./0_prep_ngsfilters/" + primernames):
        print("Primers file is not found (--primers argument is incorrect).\n"
              "File must be located in project/0_prep_ngsfilters/. You don't need provide a path, only a filename.")
        exit()
    if not os.path.isfile(f"./1_ngsfilters/{library}.ngsfilter"):
        print("ngsfilter file is not found.\n"
              "File must be located in project/1_ngsfilters/. Check create_ngsfilter.py output.")
        exit()
    primers = pd.read_csv(f"./0_prep_ngsfilters/{primernames}")
    # 1. Activate obitools3 environment
    os.system(f". {obipath}/obi3-env/bin/activate") #????
    # 2. Import reads to database
    print(f"Importing {reads1} to database")
    os.system(f"obi import --fastq-input {reads1} ./{library}/reads1")
    print(f"Importing {reads2} to database")
    os.system(f"obi import --fastq-input {reads2} ./{library}/reads2")
    # ADD a quality check for .ngsfilter file!
    # 4. Add ngsfilter to database and apply it to reads
    os.system(f"obi import --ngsfilter ./1_ngsfilters/{library}.ngsfilter ./{library}/ngsfilter")
    os.system(f"obi alignpairedend -R ./{library}/reads2 ./{library}/reads1 ./{library}/aligned_reads")
    print(f"Alignment complete. Filter out unaligned sequences")
    os.system(f"obi grep -a mode:alignment ./{library}/aligned_reads ./{library}/good_sequences")
    os.system(f"obi ngsfilter -t ./{library}/ngsfilter -u ./{library}/unidentified_sequences ./{library}/good_sequences ./{library}/identified_sequences")
    # 4. Align paired-end reads and filter out unaligned reads
    # 5. Split dataset by loci and filter it.
    print("Alignment completed and filtered. Splitting data by loci and output to .tab files")
    for loci in primers.iloc[:, 0]:
        print(f"Proceeding locus {loci}")
        os.system(f"obi grep -a experiment:{loci} ./{library}/identified_sequences  ./{library}/{loci}")  # Split data
        os.system(f"obi uniq -m sample ./{library}/{loci} ./{library}/{loci}_uniq")  # Remove PCR duplicates
        os.system(f"obi grep -p \"sequence[\'COUNT\']>=1\" ./{library}/{loci}_uniq ./{library}/{loci}_count")  # Remove all loci with number of reads < 9
        os.system(f"obi annotate --length -k COUNT -k MERGED_sample ./{library}/{loci}_count ./{library}/{loci}_cleaned") # Remove unnecessary columns
        os.system(f"obi export --tab-output  ./{library}/{loci}_cleaned > ./{library}_tab/{library}_{loci}.uniq.tab")  # library_loci.uniq.tab

if __name__ == '__main__':
    main()
