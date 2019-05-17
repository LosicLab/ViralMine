ViralMine
=============

The ViralMine pipeline is a collection of bioinformatic tools designed to process the unmapped reads from RNA-Seq or DNA-Seq alignmers and recover large viral sequence contigs of interest, as specified by the user. 

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

[TBD]

(Current version: v0.3)