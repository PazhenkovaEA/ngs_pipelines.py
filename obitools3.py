import click
import pandas as pd
import os

@click.command()
@click.option("--project", help="Path to project folder")
@click.option("--library", help="Library name")
@click.option("--reads1", help="Forward .fastq reads. Filename, should be located in Project folder")
@click.option("--reads2", help="Reverse .fastq reads. Filename, should be located in Project folder")
@click.option("--primers", help="Primers for each loci, path to .csv file")


def main(project, library, reads1, reads2, primers):
    #os.chdir(project)
    #"project" is a folder and "library" is an Illumina run name. Obi3 DMS will be named as "library"
    if not os.path.exists(f"{project}/{library}_tab"):
        os.mkdir(f"{project}/{library}_tab")
    if not os.path.isfile(primers):
        print("Primers file is not found (--primers argument is incorrect).")
        exit()
    if not os.path.isfile(f"{project}/ngsfilters/{library}.ngsfilter"):
        print("ngsfilter file is not found.\n"
              f"File must be located in {project}/ngsfilters. Check create_ngsfilter.py output.")
        exit()
    primernames = pd.read_csv(primers)
    primernames = primernames[["locus", "primerF", "primerR"]].drop_duplicates()
    # 2. Import reads to database
    print(f"Importing {reads1} to database")
    os.system(f"obi import --fastq-input {reads1} {project}/{library}/reads1")
    print(f"Importing {reads2} to database")
    os.system(f"obi import --fastq-input {reads2} {project}/{library}/reads2")
    # ADD a quality check for .ngsfilter file!

    # 3. Align paired-end reads and filter out unaligned reads
    os.system(f"obi alignpairedend -R {project}/{library}/reads2 {project}/{library}/reads1 {project}/{library}/aligned_reads")
    print(f"Alignment complete. Filter out unaligned sequences")
    os.system(f"obi grep -a mode:alignment {project}/{library}/aligned_reads {project}/{library}/good_sequences")
    # 4. Add ngsfilter to database and apply it to reads
    os.system(f"obi import --ngsfilter {project}/ngsfilters/{library}.ngsfilter {project}/{library}/ngsfilter")
    os.system(f"obi ngsfilter -t {project}/{library}/ngsfilter -u  {project}/{library}/unidentified_sequences {project}/{library}/good_sequences {project}/{library}/identified_sequences")

    # 5. Split dataset by loci and filter it.
    print("Alignment completed and filtered. Splitting data by loci and output to .tab files")
    for loci in primernames.iloc[:, 0]:
        print(f"Proceeding locus {loci}")
        os.system(f"obi grep -a experiment:{loci} {project}/{library}/identified_sequences  {project}/{library}/{loci}")  # Split data
        os.system(f"obi uniq -m sample {project}/{library}/{loci} {project}/{library}/{loci}_uniq")  # Remove PCR duplicates
        #os.system(f"obi grep -p \"sequence[\'COUNT\']>=1\" {project}/{library}/{loci}_uniq {project}/{library}/{loci}_count")  # Remove all loci with number of reads < 9
        os.system(f"obi annotate --length -k COUNT -k MERGED_sample {project}/{library}/{loci}_uniq {project}/{library}/{loci}_cleaned") # Remove unnecessary columns
        os.system(f"obi export --tab-output  {project}/{library}/{loci}_cleaned > {project}/{library}_tab/{library}__{loci}.uniq.tab")  # library__loci.uniq.tab

if __name__ == '__main__':
    main()
