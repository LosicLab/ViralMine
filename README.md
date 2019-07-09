ViralMine
=============

The ViralMine pipeline is a collection of bioinformatic tools designed to process the unmapped reads from RNA-Seq or DNA-Seq aligners and recover large viral sequence contigs of interest, as specified by the user. 

As the original purpose of the program was to recover contigs matching the Hepatitis B virus, additional functionality for genotyping found viral contigs has been included, using a similar method to the NCBI virus genotyping web tool: https://www.ncbi.nlm.nih.gov/projects/genotyping/formpagex.cgi

To reduce dependencies and program complexity, the software currently requires the user to have already aligned their sequencing reads using STAR, HISAT2, bwa, etc. *Please be sure to adjust your aligner parameters to keep Unmapped reads in distinct output(s)!*

Please see the docs for more complete information: [URL TO BE INCLUDED]


## Downloads and Prerequisites ##

To use ViralMine, download a clone of this repositiory (git clone https://github.com/LosicLab/ViralMine.git) and add the `ViralMine` folder to your PATH (you can do this in linux/MacOS by modifying your .bashrc file, with the command `export PATH=$PATH:path/to/ViralMine`)

Additionally, you will need to have the following programs installed and added to your path before you can mine your first sample:

```
BLAST (recommend latest version)
python (require v3.6.x or later)
cd-hit (recommend v4.6.x or later)
TRINITY (see here: https://github.com/trinityrnaseq/trinityrnaseq/wiki)
```

**Note**: because of the large memory requirements of TRINITY's inchworm assembler, it is *highly* recommended that you run ViralMine on a computing cluster with at least 40GB of RAM. If this is not possible, you will need to adjust the default JM memory parameters in the inchworm step of the `ViralMine.sh` script.   


## Quick Start Guide ##

Run parameters are located in the first code block of the `ViralMine.sh` file, and should be modified to point to the required alignment directory, viral sequence fasta, viral BLAST database output directory, and more. Detailed parameter descriptions can be found in the docs (See above).

Files Required:

```
1. Viral reference sequences in fasta format (or, can use the provided HBVdb fasta or BLASTdb [see 'HBV_Ref_dbs'])
2. Unmapped reads in fastq format (eg. for single end STAR alignments, these will be 'Unmapped.reads.out'. Check your aligner manual for more information)
```

Once you have set your parameters in the script file, execute the pipeline (for single end reads) using:

```
ViralMine.sh Unmapped.reads.out
```

Replacing `Unmapped.reads.out` with the file name of the unmapped reads file depending on your aligner.

For paired end sequencing, just specify the additional mate pair:

```
ViralMine.sh Unmapped.reads.out1 Unmapped.reads.out2
```

Please see the docs for additional information.


## Output Files ##

Each step of the pipeline will produce several output files, and depending on the size of your unmapped read fastqs, you should expect to use 5-15GB of storage. Key output files are summarized below:

1. `viral_matched_contigs.fa`: A fasta file containing all the inchworm contigs that matched viral reference sequences
2. `viral_alignment.tsv`: The BLAST output with scores of which contigs matched which viral sequences. This can be used to identify which contigs matched to which viral species/viral reference.
3. `[sample_name]_scores.txt`: Will only be generated if the HBV genotyping flag has been selected. This will contain the bitscores by genotype for the BLAST window alignment, and can be used to genotype the HBV of a patient, or characterize a mixed genotype.
4. `[sample_id]_viral_GT.tsv`: Will only be generated if the HBV genotyping flag has been selected. This contains the calculated dominant viral genotype given the BLAST results. 

## Known Issues ##

1. Errors will be returned when there are no viable contigs produced by inchworm or after filtering, as no further steps can be completed. Additional error will be thrown if directory structures are not maintained -- for this reason *it is important to use absolute paths where possible in the pipeline parameters!*

2. Viral genotyping of HBV is tricky, and is dependent upon the viral reference database you are using. Current runs with the software have shown that Refseq and HBVdb sequences are similar enough that overall viral genotype should be called consistently, but bitscore/mixed genotype calls will differ depending on reference source. 


(Current version: v0.3)
