# ChimeraKiller
In-development script for identifying chimeric transcripts in a de novo assembled transcriptome


Suggestion for creating conda environment (for version 0.7.3 and up which use gatk4 and the pathos multiprocessing library):


conda create -n chimerakiller_env python=3.6 biopython bwa samtools bedtools picard pandas matplotlib scipy pysam gatk4 pathos openssl=1.0 minimap2

conda activate chimerakiller_env



## What does ChimeraKiller do?
ChimeraKiller was written with the goal of providing a quality control and filtering step for de novo assembled transcriptomes generated from short read illumina data. ChimeraKiler primarily works by assessing read mapping of the reference's raw reads against the transcript set to attempt to identify likely misassembled and/or chimeric transcripts. It then filters these sequences into 'good', 'bad', and 'bad_low' folders based on failed filters and read coverage. ChimeraKiller will also output PDF coverage plots for each sequence which can be viewed by the user to observationally evaluate and correct ChimeraKiller's sorting if they wish.

ChimerKillers procedure is as follows:

In the first step, reads are mapped to the transcripts with BWA or, more recently, minimap2. After read mapping, reads are filtered with samtoools and gatk such that only those with an allowable number of mismatches are kept (the default 0 mismatches). Coverage is then calculated weith bedtools for every site in the set of transcripts. 

The first filtering step is based on coverage criteria. Any transcripts with sites that have 0 coverage are identified as 'bad' and sorted accordingly. If the defaults are used, this will identify any transcript without an exact read match across its full length as 'bad'. The logic here being that the sequence cannot be reliable if it cannot be unambiguously assembled by the available reads.

Sequences that pass the first filtering step continue on to site-based filtering of unevenly distributed reads. This filtering step aims to identify transcripts where read mappings do not span uniformally across all sites. For a given transcript where reads map randomly across the sequence length, we would expect the matching site to occur at similarly random positions in reads (Figure 1).
