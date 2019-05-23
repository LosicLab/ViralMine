ViralMine: Viral Sequence Miner Tool
====================================

| Created in the Losic Lab at Ichan Institute of Genetics and Genomics, 
| Mount Sinai School of Medicine, New York, NY

Introduction
------------

The **ViralMine** pipeline is a collection of bioinformatic tools designed to process the unmapped reads from RNA-Seq or DNA-Seq aligners and recover large viral sequence contigs of interest, as specified by the user. 

To reduce dependencies and program complexity, this software currently requires the user to have already aligned their sequencing reads using STAR, HISAT2, bwa, etc. *Please be sure to adjust your aligner parameters to keep Unmapped reads in distinct output(s)!*. Please see your aligner-of-choice's manual for how to do this.

As the original purpose of the program was to recover contigs matching the Hepatitis B virus, additional functionality for genotyping found viral contigs has been included, using a similar method to the NCBI virus genotyping web tool: https://www.ncbi.nlm.nih.gov/projects/genotyping/formpagex.cgi


Installation
------------

**ViralMine** is designed to be run from the MacOS command line or linux shell environment, an uses a combination of Bash and python. 

Because of the large RAM requirements of component programs for this pipeline (TRINITY's Inchworm) it is *highly* recommended that *ViralMine* be installed on a compute cluster or a machine with at least 30-40GB of RAM available. If this is not possible for you, you will need to adjust the trinity jelly-fish memory parameters to a memory size your machine can allocate ("-JM [Memory]").

To install and set up **ViralMine**, from the command line:

1. Navigate to the directory where you want to install the pipeline
2. Clone *ViralMine* into your local repositiory: 
	
		``git clone https://github.com/LosicLab/ViralMine.git``

3. Make sure the bash and python scripts have execute permissions:
	
		``chmod -x ViralMine.sh``
		``chmod -x scripts/*.py``

4. Add the cloned repository to your PATH:
	
		``export PATH=$PATH:{path/to/ViralMine}``

At this point, check that **ViralMine** has been added to your path using `echo $PATH`.

Additionally, you will need to have the following programs installed and added to your path before you can mine your first sample:

- BLAST >= 2.6.2 (recommend latest version)
- python >=3.6.x (recommend latest version)
- cd-hit >=4.6.x (recommend latest version)
- TRINITY (https://github.com/trinityrnaseq/trinityrnaseq/wiki)


Running ViralMine
=================


Setup: BLAST db
---------------

**ViralMine** requires the user to have generated a BLAST nucleotide database of viral reference sequences before samples can be processed. This can be done by using a fasta file of collected sequences and the command line BLAST tools' ``makeblastdb`` command (see the BLAST manual here for more information: https://www.ncbi.nlm.nih.gov/books/NBK279690/). Make sure the directory for you compiled database is in a location **ViralMine** can read from.

If you are mining HBV viral sequences, a nucleotide reference database has already been compiled for HBV genotypes `A, B, C, D, F,` and `G` from the HBVdb reference (https://hbvdb.ibcp.fr/HBVdb/HBVdbIndex) and included in the **ViralMine** download (see the ``HBV_Ref_dbs/`` directory).


Setup: Alignment
----------------

Prior to running **ViralMine**, users will need to generate the aligned and unaligned read files for samples they wish to process. There is no specific aligner requirement, but STAR has been the validated aligner of choice for RNA-Seq reads, and BWA for DNA-Seq reads. Check your aligner of choice's manual to determine what parameter flag must be set to generate *unmapped reads*.

**ViralMine** does not require any files other than the Unmapped Reads file(s) from the aligner output. 


Setup: ViralMine parameter configuration:
-----------------------------------------

**ViralMine** is currently configured to run on one single or paired-end RNA/DNA sample at a time**. Parameters are specified in the opening block of the ``ViralMine.sh`` script, seen below:
	
	```
	# Parameters:
	Dir="path/to/unmapped_reads.out/" #Path to the location of your alignment output files
	seq_type="paired" # Select "paired" or "single" end sequencing (to select how many fastqs to expect)
	Viral_Genome="path/to/input.genome.fa" # Input fasta containing viral reference sequence(s)
	viral_db="path/to/viral/blastn_db/viral.db" # Where blast database for viral reference sequences will be output (can be substituted with reference db included in git repo)
	sample_id="sample_name" # Name you want to give the sample (JobID)
	contig_size_filter=200 #length flag below which putative viral contigs will be removed  
	gt_virus=1 # 0 or 1 (No or Yes) to specify if viral contigs should be genotyped (built for HBV, currently)
	```

A full list of the parameters and their options is discussed in the table below (see **Parameter Explanations**).


It is recommended that for *each sample that you run ViralMine with, you make a uniquely named copy of the ``ViralMine.sh`` executable. For example, if you have samples ``Sample_1`` and ``Sample_2``, you would create ``ViralMine_S1.sh`` and ``ViralMine_S2.sh`` in the ``~/ViralMine`` directory. You do **NOT** need to copy any of the files in ``scripts/``.


Running ViralMine
-----------------

Once you have generated your unmapped read files, built your viral reference database, and specified your input files and reference directories in the ``ViralMine.sh`` parameters, you can execute the script for single-end samples using:

```
ViralMine.sh Unmapped.Reads.out
```

or for paired-end samples using

```
ViralMine.sh Unmapped.Reads.out1 Unmapped.Reads.out2
```

Script processes and error messages are configured to be printed to standard out. If you would like to capture these in a log file, run:


```
ViralMine.sh Unmapped.Reads.out > {path/to/log_dir}/log.out
```


Parameter Explanations 
======================

[TBD]



HELP
====

If you have further questions, you can email me at adrian.bubie@mssm.edu
