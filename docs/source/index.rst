ViralMine: Viral Sequence Miner Tool
====================================

| Created in the Losic Lab at Ichan Institute of Genetics and Genomics, 
| Mount Sinai School of Medicine, New York, NY

Introduction
------------

The **ViralMine** pipeline is a collection of bioinformatic tools designed to process the unmapped reads from RNA-Seq or DNA-Seq aligners and recover large viral sequence contigs of interest, as specified by the user. 

To reduce dependencies and program complexity, this software currently requires the user to have already aligned their sequencing reads using STAR, HISAT2, bwa, etc. *Please be sure to adjust your aligner parameters to keep Unmapped reads in distinct output(s)!* Please see your aligner-of-choice's manual for how to do this.

As the original purpose of the program was to recover contigs matching the Hepatitis B virus and Human Papillomavirus, additional functionality for genotyping found viral contigs has been included, using a similar method to the NCBI virus genotyping web tool: https://www.ncbi.nlm.nih.gov/projects/genotyping/formpagex.cgi

Additionally, **ViralMine** can generate putative single-virus coinfections for HBV and HPV infected samples, and can produce viral gene expression counts (supporting reads) for recovered viral sequences without additional steps.


Installation
------------

**ViralMine** is designed to be run from the command line environment, an uses a combination of Bash and python. 

Because of the large RAM requirements of component programs for this pipeline (TRINITY's Inchworm) it is *highly* recommended that *ViralMine* be installed on a compute cluster or a machine with at least 30-40GB of RAM available. If this is not possible for you, you will need to adjust the trinity maximum memory parameters to a memory size your machine can allocate (``--max_memory [Gb]``).

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
- TRINITY (require 2.8.x or later; https://github.com/trinityrnaseq/trinityrnaseq/wiki)

Finally, *before running the pipeline the first time*, you will need to adjust the ``vm_loc`` parameter in the ``parameters.txt`` file to the location where the **ViralMine** git repository has been cloned you your machine. This will allow the script to find the appropriate files required to run the pipeline.

Running ViralMine
=================

Setup: Alignment
----------------

Prior to running **ViralMine**, users will need to generate the aligned and unaligned read files for samples they wish to process. There is no specific aligner requirement, but STAR has been the validated aligner of choice for RNA-Seq reads, and BWA for DNA-Seq reads. Check your aligner of choice's manual to determine what parameter flag must be set to generate *unmapped reads*.

**ViralMine** does not require any sequencing files other than the Unmapped Reads file(s) from the aligner output. 


Setup: BLAST db
---------------

**ViralMine** requires the user to generate a BLAST nucleotide database of viral reference sequences before samples can be processed. This can be done by using a fasta file of collected sequences and the command line BLAST tools' ``makeblastdb`` command (see the BLAST manual here for more information: https://www.ncbi.nlm.nih.gov/books/NBK279690/). This is handled in the first step of the pipeline, but if the user already has a BLAST database compiled, they can skip this step using the ``Exisiting_Blastdb`` flag.

If you are mining HBV or HPV viral sequences, a nucleotide reference database has already been compiled for HBV genotypes ``A, B, C, D, F,`` and ``G`` from the HBVdb reference (https://hbvdb.ibcp.fr/HBVdb/HBVdbIndex) and included in the **ViralMine** download (see the ``HBV_Ref_dbs/`` directory). The script defualt is to pointed to this location.

If you are mining HPV viral reads, a nucleotide reference database is also available (``HPV_Ref_dbs/``) that includes the 15 oncogenic genotypes associated with the vast majority of cervical cancers (*Castellsague, X. Gynecol Oncol. 2008*). Genome reference sequences for these genotypes were compiled from the most recent published sequence by RefSeq.


Setup: ViralMine parameter configuration:
-----------------------------------------

**ViralMine is currently configured to run on one single or paired-end RNA/DNA sample at a time**. Parameters are specified in the ``parameters.txt`` file, that must be passed as the first variable when running the script. The parameters are listed below:

::

	Dir="path/to/unmapped_reads.out/"
	sample_id="sample_name"
	seq_type="paired" 
	Exisiting_Blastdb="Yes"
	Viral_Genome="path/to/input.genome.fa" 
	viral_db="/path/to/ViralMine/HBV_Ref_dbs/HBVdb/HBVdb_all_gt"
	contig_size_filter=100  
	threshold=0.1
	gt_virus="hbv"
	gene_exp="No"
	viral_gene_db="/path/to/ViralMine/HBV_Ref_dbs/HBV_gene_db/GenesHBV"
	 
	

A full list of the parameters and their options is discussed in the table below (see **Parameter Explanations**).


It is recommended that for *each sample* that you run **ViralMine** with, you make a uniquely named copy of the ``parameters.txt`` file. For example, if you have samples ``Sample_1`` and ``Sample_5``, you would create ``parameters_S1.txt`` and ``parameters_S5.txt`` in your working directory. You do **NOT** need to copy any of the files in ``scripts/``.


Running ViralMine
-----------------

Once you have generated your unmapped read files, built your viral reference database, and specified your input files and reference directories in the ``parameters.txt`` file, you can execute the script for single-end samples using:

``ViralMine.sh parameters.txt [Unmapped.Reads.out] > VM_log.out``

or for paired-end samples using

``ViralMine.sh parameters.txt [Unmapped.Reads.out1] [Unmapped.Reads.out2] > VM_log.out``

Script processes and error messages are configured to be printed to standard out, which will be captured in the specified log file. For verbose output to standard out, execute the command without the log file specification:

``ViralMine.sh parameters.txt [Unmapped.Reads.out]``




Parameter Explanations 
======================

:``Dir``:	Directory location of the unmapped reads file(s), as well as the location where the output files will be published. It is highly recommended that the absolute path be used.
:``sample_id``:	Sample name or sample ID. This will be used to name the outfiles. Please note that the "_" character is restricted and must not be used.
:``seq_type``:	Either "paired" (default) or "single". Flag used to specify whether paired or single-ended sequencing was used, and to specify how many unmapped reads files the script should expect.
:``Exisiting_Blastdb``:	Either "No" or "Yes". This indicates whether or not a new nucleotide BLASTdb needs to be built from the passed in viral reference fasta. Default is Yes, as the default HBV Reference database has been included in the ViralMine download.
:``Viral_Genome``:	Filepath to fasta containing viral reference sequences to build a new nucleotide BLAST database. Will be ignored it "Exisiting_Blastdb" is 1.
:``viral_db``:	Path to either existing viral reference BLASTdb OR path to where the viral database should be created.
:``contig_size_filter``:	Integer value, specifying the smallest contig size to keep when aligning agains the viral references. Default size is 100bp (what we have found to work well for HBV).
:``gt_virus``:	Flag for running contig genotyping, either "hbv","hpv", or "none". **Note** this should ONLY be run if you are trying to genotype HBV or HPV contigs with the included reference databases! This will most likely fail or provide useless results for other viruses. Default is "hbv", as the expected default reference database is HBV. 
:``threshold``:	Fraction of the total viral bitscore for a patient that must be exceeded for a given genotype of HBV/HPV to be included in the list of coinfections in the coinfection output file. The default is 0.1 or 10% of the bitscore. Will only be used if Genotyping is attempted.
:``gene_exp``:	Flag to indicate whether or not to generate viral gene level expression count matrix for HPV or HBV infected samples. Either "No" or "Yes" (default is "No")
:``viral_gene_db``:	Path to the HBV or HPV viral gene nucleotide BLAST reference database. *NOTE*: this MUST be either the HPV or HBV reference database provided in the ViralMine download to function correctly. HBV reference db is "HBV_Ref_dbs/HBV_gene_db/GenesHBV" and HPV is "HPV_Ref_dbs/HPV_gene_db/genesHPV".




Output Files
============

Output files will be saved in the same directory as the input directory where your unmapped reads files are located. Each step of the pipeline will produce several output files, and depending on the size of your unmapped read fastqs, you should expect to use 5-10GB of storage. Key output files are summarized below:

1. ``viral_matched_contigs.fa``: A fasta file containing all the inchworm contigs that matched viral reference sequences
2. ``viral_alignment.tsv``: The BLAST output with scores of which contigs matched which viral sequences. This can be used to identify which contigs matched to which viral species/viral reference.
3. ``[sample_id]_scores.txt``: Will only be generated if the HBV or HPV genotyping flag has been selected. This will contain the bitscores by genotype for the BLAST window alignment, and can be used to genotype tumor viral infection of a patient, or characterize a mixed genotype.
4. ``[sample_id]_viral_GT.tsv``: Will only be generated if the genotyping flag is on for HPV or HBV. This will contain the calculated dominant genotype for the tumor across patient contigs by summation of all genotype specific bitscores from across the window BLAST.
5. ``[sample_id]_viral_Coinf_GT.tsv``: Will only be generated if the HBV/HPV genotyping flag has been selected. This will contain a comma separated list of the viral genotypes that the patient's tumor is coinfected with. It is NOT ordered. If only one genotype is listed for a patient, this indicates that only one genotype passed the coinfection threshold.
6. ``[sample_id]_ReadsPerViralGene.tab``: Will only be generated if the `gene_exp flag` is "Yes", and the HBV/HPV genotyping flag has been selected. This will calculate, based on the recovered viral contigs, the number of reads supporting HBV/HPV genes for the patient.



HELP
====

If you have further questions, you can email me at adrian.bubie@mssm.edu
