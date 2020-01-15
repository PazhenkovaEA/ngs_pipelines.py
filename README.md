# Analysis of NGS data for large-scale conservation genetics projects

This repo contains python scripts for demultiplexing NGS data, allele calling and prepares data for analysis with The Database. Main steps of demultiplexing are based on OBITools3: https://git.metabarcoding.org/obitools/obitools3/blob/master/README.md
## Experiment design

Samples are stored in 96-wells **aliquot plates**.
The aliquot plate can be pipetted to one or several (if we need quality control replicates) **primer plates**.
Primer plates (or, better to say, Tag plates?) contain unique combination of tags (e.g. acacacac:acacacac) for each position.
This tags are appended to the all forward and reverse primers.
We perform a multiplex PCR, in each well we got PCR product for all loci of interest, marked with the position-specific tag combination.
After sequencing, we need to detect all the samples and loci and then perform allele calling.
## Installation

OS X & Linux
Install OBITools3 following to [developers’ instruction](https://git.metabarcoding.org/obitools/obitools3/wikis/Installing-the-OBITools3)

python3, python3-venv, git, Cmake need to be installed. Also python packages Click==7.0, numpy==1.17.4, pandas==0.25.3 should be installed to obi3 environment. You can install it manually or use the following command.
```sh
git clone https://github.com/PazhenkovaEA/ngs_pipelines.py
cd ngs_pipelines.py
pip install -r ngs_req.txt
```

## Usage example
### Commands 
1. Create ngsfilter file 
```sh
python3 create_ngsfilter.py --project=./UA1
--plates=./UA1/AliquotP.xlsx --tags=./UA1/UA_tagscombo.csv –-primernames=./UA1/UA_primers.csv --aliquotplates=AP
```
Produces files ./UA1/ngsfilters/UA1.ngsfilter (for searching of primers and tags with OBITools) and ./UA1/results/UA1_sample_positions.txt (used for quality control).

2. Obitools part

Don’t forget to activate obi3 environment
```sh
. ~/obitools3/obi3-env/bin/activate
```
If your project contains several library, you should run this script separately for each library.
```sh
python obitools3.py --project=./UA1 --library=UA1 --reads1=UA_F.fastq --reads2=UA_R.fastq –primers=./UA1/UA_primers.csv
```

Creates a obitools3 database and a directory library_tab with found sequences for each loci.

3. Allele calling

```sh
python callAllele.py --project=./UA1 --primers=./UA1/UA_primers.csv 
```

### Input description
* Illumina paired-end reads in .fastq format
Examples of the input files can be found in /0_prep_ngsfilters directory.
Copy this directory to yours project folder and use input files as a templates.
In the example we have the library UA1, Aliquot Plate A1 with samples bear1 - bear12 and Aliquot Plate A2 with samples bear13 - bear24, Primer Plates PP1 and PP2, microsatellite UA_03, sex marker ZF (snp) and mitochondrial marker COI.
*All tags, marker and primer sequences in the example are not real and usable for educational purposes only.*
* Library, aliquot plates and primer plates mapping (keep the column names, table in .xlsx)
Even if you don't use Aliquot Plates and Primer Plates, this file is obligate.
You can pick any Aliquot Plate name.

 **Please, avoid . and / symbols in loci and sample names**.

Example:

Library_BC| Aliquot Plate | Primer Plate
--- | --- | ---
UA1|  A1 | PP1
UA1|  A2 | PP2
* Samples localization on plates (separate table in .xls(x) for each aliquot plate).
The most important file(s).
Files must be named **AP_library_AliquotePlate** (e.g., AP_UA1_A1.xls, AP_UA1_A2.xls) and located in one folder.
Three columns are obligate: "SPositionBC" with the sample names, "TPositionId" with sample location in format A1, B1, C1...H12 and "SPositionId" with sample location in numeric format 1, 2..96.
If you don't use Aliquot Plates, you will need the only one file.
Example:

SPositionBC| TPositionId | SPositionId
--- | --- | ---
bear1|   A1| 1
bear1|  B1  | 2
...|  ... | ...
bear12|  H12 | 96
* Combination of tags in .csv format (columns: position,primer plate 1,... primer plateN).
If you have the only one primer plate, call the column PP1 anyway.
Example:

position| PP1 | PP2
--- | --- | ---
1|  acacacac:acacacac |caggctaa:tgagccta
2|  acacacac:acagcaca |caggctaa:tgagcctt
..|  ...|...
96|  gactgatg:gatcgcga |gaggacta:tcagtcga
* List of primers, motifs and reference sequences for each locus in **.csv** format (locus, primerF, primerR, type, motif, sequence).
loci names have to contain letters.
Allowed types: microsat, snp and mt.
In case if the marker is biallelic or microsatellite contains different motifs, you need to provide all sequences or motifs with the same locus name.
Sequence (in the "sequence" column) should be provided without primers.

Example:

locus | primerF | primerR | type | motif |sequence
--- | --- | --- | --- | --- | ---
UA_03|  gctcccataac |gctcccataac | microsat | acac |
ZF|  cataacgctcc |taacgctcccataac | snp | | agag........tatac |
ZF|  cataacgctcc |taacgctcccataac | snp | | agag........tagac |
COI|  gaatcgccacc |acatcaag | mt | | CT.....TAAACTATTCCCTG |
### OBITools pipeline description
Here is a pseudocode for Obi3 pipeline with some comments.
1. Import reads to the database, align it and filter out unaligned reads

```sh
obi import --fastq-input reads1 /library/reads1
obi import --fastq-input reads2 /library/reads2
obi alignpairedend -R /library/reads2 /library/reads1 /library/aligned_reads
obi grep -a mode:alignment /library/aligned_reads /library/good_sequences
```
2. Import ngsfilter and apply it to the aligned reads
```sh
obi import --ngsfilter file.ngsfilter /library/ngsfilter
obi ngsfilter -t /library/ngsfilter -u /library/unidentified_sequences /library/good_sequences /library/identified_sequences

```
3. 
Split data by loci, remove duplicates (add COUNT column), remove unnecessary columns and export data to .tab file.
```sh
for loci in primers.csv;
do
obi grep -a experiment:$loci /library/identified_sequences  /library/$loci;
obi uniq -m sample /library/$loci /library/$loci_uniq;
obi annotate --length -k COUNT -k MERGED_sample /library/$loci_uniq /library/$loci_cleaned;
obi export --tab-output  /library/$loci_cleaned > library__$loci.uniq.tab;
done
```
## References

The previous version of this pipeline was developed by [Roman Luštrik](https://github.com/romunov/ngs_pipelines)
The pipeline follows a methodology described in:
De Barba, M., Miquel, C., Lobréaux, S., Quenette, P.Y., Swenson, J.E. and Taberlet, P. (2017), High‐throughput microsatellite genotyping in ecology: improved accuracy, efficiency, standardization and success with low‐quantity and degraded DNA. Mol Ecol Resour, 17: 492-507. doi:10.1111/1755-0998.12594
