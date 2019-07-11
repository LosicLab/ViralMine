ViralMine: Viral Sequence Miner Tool
====================================

| Created in the Losic Lab at Ichan Institute of Genetics and Genomics, 
| Mount Sinai School of Medicine, New York, NY

Introduction
------------

The **ViralMine** pipeline is a collection of bioinformatic tools designed to process the unmapped reads from RNA-Seq or DNA-Seq aligners and recover large viral sequence contigs of interest, as specified by the user. 

To reduce dependencies and program complexity, this software currently requires the user to have already aligned their sequencing reads using STAR, HISAT2, bwa, etc. *Please be sure to adjust your aligner parameters to keep Unmapped reads in distinct output(s)!* Please see your aligner-of-choice's manual for how to do this.

As the original purpose of the program was to recover contigs matching the Hepatitis B virus and Human Papillomavirus, additional functionality for genotyping found viral contigs has been included, using a similar method to the NCBI virus genotyping web tool: https://www.ncbi.nlm.nih.gov/projects/genotyping/formpagex.cgi


Installation
------------

**ViralMine** is designed to be run from the MacOS command line or linux shell environment, an uses a combination of Bash and python. 

Because of the large RAM requirements of component programs for this pipeline (TRINITY's Inchworm) it is *highly* recommended that *ViralMine* be installed on a compute cluster or a machine with at least 30-40GB of RAM available. If this is not possible for you, you will need to adjust the trinity jelly-fish memory parameters to a memory size your machine can allocate (``-JM [Memory]``).

To install and set up **ViralMine**, from the command line:

1. Navigate to the directory where you want to install the pipeline
2. Clone *ViralMine* into your local repositiory: 
	
		``git clone https://github.com/LosicLab/ViralMine.git``

3. Make sure the bash and python scripts have execute permissions:
	
		``chmod -x ViralMine.sh scripts/*.py``

4. Add the cloned repository to your PATH:
	
		``export PATH=$PATH:{path/to/ViralMine}``

At this point, check that **ViralMine** has been added to your path using ``echo $PATH``. **NOTE:** this will only add the path temporarily (until you log out of the terminal); to permanantly add **ViralMine** to your path, add the command in 4. to your '.bashrc' file, and run ``source ~/.bashrc``. 

Additionally, you will need to have the following programs installed and added to your path before you can mine your first sample:

- BLAST >= 2.6.2 (recommend latest version)
- python >=3.6.x (recommend latest version)
- cd-hit >=4.6.x (recommend latest version)
- TRINITY (https://github.com/trinityrnaseq/trinityrnaseq/wiki)


Running ViralMine
=================

Setup: Alignment
----------------

Prior to running **ViralMine**, users will need to generate the aligned and unaligned read files for samples they wish to process. There is no specific aligner requirement, but STAR has been the validated aligner of choice for RNA-Seq reads, and BWA for DNA-Seq reads. Check your aligner of choice's manual to determine what parameter flag must be set to generate *unmapped reads*.

**ViralMine** does not require any files other than the Unmapped Reads file(s) from the aligner output. 



Setup: BLAST db
---------------

**ViralMine** requires the user to generate a BLAST nucleotide database of viral reference sequences before samples can be processed. This can be done by using a fasta file of collected sequences and the command line BLAST tools' ``makeblastdb`` command (see the BLAST manual here for more information: https://www.ncbi.nlm.nih.gov/books/NBK279690/). This is handled in the first step of the pipeline, but if the user already has a BLAST database compiled, they can skip this step using the ``Exisiting_Blastdb`` flag.


If you are mining HBV viral sequences, a nucleotide reference database has already been compiled for HBV genotypes `A, B, C, D, F,` and `G` from the HBVdb reference (https://hbvdb.ibcp.fr/HBVdb/HBVdbIndex) and included in the **ViralMine** download (see the ``HBV_Ref_dbs/`` directory). The script defualt is to pointed to this location.

If you are mining HPV viral reads, a nucleotide reference database is also available (``HPV_Ref_dbs/``) that includes the 15 oncogenic genotypes associated with the vast majority of cervical cancers (*Castellsague, X. Gynecol Oncol. 2008*). Genome reference sequences for these genotypes were compiled from the most recent published sequence by RefSeq.


Setup: ViralMine parameter configuration:
-----------------------------------------

**ViralMine is currently configured to run on one single or paired-end RNA/DNA sample at a time**. Parameters are specified in the opening block of the ``ViralMine.sh`` script, seen below:

::

	Dir="path/to/unmapped_reads.out/" 
	seq_type="paired" 
	Exisiting_Blastdb=1
	Viral_Genome="path/to/input.genome.fa" 
	viral_db="path/to/viral/blastn_db/viral.db"
	contig_size_filter=200  
	gt_virus="hbv"
	sample_id="sample_name" 
	

A full list of the parameters and their options is discussed in the table below (see **Parameter Explanations**).


It is recommended that for *each sample* that you run **ViralMine** with, you make a uniquely named copy of the ``ViralMine.sh`` executable. For example, if you have samples ``Sample_1`` and ``Sample_2``, you would create ``ViralMine_S1.sh`` and ``ViralMine_S2.sh`` in the ``~/ViralMine`` directory. You do **NOT** need to copy any of the files in ``scripts/``.


Running ViralMine
-----------------

Once you have generated your unmapped read files, built your viral reference database, and specified your input files and reference directories in the ``ViralMine.sh`` parameters, you can execute the script for single-end samples using:

``ViralMine.sh Unmapped.Reads.out``

or for paired-end samples using

``ViralMine.sh Unmapped.Reads.out1 Unmapped.Reads.out2``

Script processes and error messages are configured to be printed to standard out. If you would like to capture these in a log file, run:

``ViralMine.sh Unmapped.Reads.out > {path/to/log_dir}/log.out``


Parameter Explanations 
======================

:``Dir``:	Directory location of the unmapped reads file(s), as well as the location where the output files will be published. It is highly recommended that the absolute path be used.
:``seq_type``:	Either "paired" (default) or "single". Flag used to specify whether paired or single-ended sequencing was used, and to specify how many unmapped reads files the script should expect.
:``Exisiting_Blastdb``:	Either "No" or "Yes". This indicates whether or not a new nucleotide BLASTdb needs to be built from the passed in viral reference fasta. Default is Yes, as the default HBV Reference database has been included in the ViralMine download.
:``Viral_Genome``:	Filepath to fasta containing viral reference sequences to build a new nucleotide BLAST database. Will be ignored it "Exisiting_Blastdb" is 1.
:``viral_db``:	Path to either existing viral reference BLASTdb OR path to where the viral database should be created.
:``contig_size_filter``:	Integer value, specifying the smallest contig size to keep when aligning agains the viral references. Default size is 100bp (what we have found to work well for HBV).
:``gt_virus``:	Flag for running contig genotyping, either "hbv","hpv", or "none". **Note** this should ONLY be run if you are trying to genotype HBV or HPV contigs with the included reference databases! This will most likely fail or provide useless results for other viruses. Default is "hbv", as the expected default reference database is HBV. 
:``sample_id``:	Sample name or sample ID. This will be used to name the genotyping outfile.


Output Files
============

Each step of the pipeline will produce several output files, and depending on the size of your unmapped read fastqs, you should expect to use 5-15GB of storage. Key output files are summarized below:

1. ``viral_matched_contigs.fa``: A fasta file containing all the inchworm contigs that matched viral reference sequences
2. ``viral_alignment.tsv``: The BLAST output with scores of which contigs matched which viral sequences. This can be used to identify which contigs matched to which viral species/viral reference.
3. ``[sample_id]_scores.txt``: Will only be generated if the HBV or HPV genotyping flag has been selected. This will contain the bitscores by genotype for the BLAST window alignment, and can be used to genotype tumor viral infection of a patient, or characterize a mixed genotype.
4. ``[sample_id]_viral_GT.tsv``: Will only be generated if the genotyping flag is on for HPV or HBV. This will contain the calculated dominant genotype for the tumor across patient contigs by summation of all genotype specific bitscores from across the window BLAST.


HELP
====

If you have further questions, you can email me at adrian.bubie@mssm.edu
