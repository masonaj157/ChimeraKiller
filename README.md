# ChimeraKiller
In-development script for identifying chimeric transcripts in a de novo assembled transcriptome


Suggestion for creating conda environment (for version 0.7.3 and up which use gatk4 and the pathos multiprocessing library):


conda create -n chimerakiller_env python=3.6 biopython bwa samtools bedtools picard pandas matplotlib scipy pysam gatk4 pathos openssl=1.0 minimap2

conda activate chimerakiller_env


## What does ChimeraKiller do?
ChimeraKiller was written with the goal of providing a quality control and filtering step for de novo assembled transcriptomes generated from short read illumina data. ChimeraKiler primarily works by assessing read mapping of the reference's raw reads against the transcript set to attempt to identify likely misassembled and/or chimeric transcripts. It then filters these sequences into 'good', 'bad', and 'bad_low' folders based on failed filters and read coverage. ChimeraKiller will also output PDF coverage plots for each sequence which can be viewed by the user to observationally evaluate and correct ChimeraKiller's sorting if they wish.

ChimeraKiller's protocol is as follows:

In the first step, reads are mapped to the transcripts with BWA or, more recently, minimap2. After read mapping, reads are filtered with samtoools and gatk such that only those with an allowable number of mismatches are kept (the default 0 mismatches). Coverage is then calculated weith bedtools for every site in the set of transcripts. 

The first filtering step is based on coverage criteria. Any transcripts with sites that have 0 coverage are identified as 'bad' and sorted accordingly. If the defaults are used, this will identify any transcript without an exact read match across its full length as 'bad'. The logic here being that the sequence cannot be reliable if it cannot be unambiguously assembled by the available reads.

Sequences that pass the first filtering step continue on to site-based filtering of unevenly distributed reads. This filtering step aims to identify transcripts where read mappings do not span uniformally across all sites. For a given transcript where reads map randomly across the sequence length, we would expect the a given site to occur at similarly random positions in mapped reads. In contrast, if a sequence is chimeric, we would expect very few reads to 'span' a given site and the site of interest would disproportionaly occur at the ends of reads. We can identify these instances by checking whether the number of bases on either side of the test site in a read shows a strong disequalibrium.

Specifically, in the the second filtering step ChimeraKiller examines each site in a sequence and calculates the average difference (base disequilibrium) in sequence length on the left and right sides of the test site in each mapped read. If reads map uniformly across sites, then the base disequilibirium is typically low. If the base disequilibrium for any site in a sequence is above 75% (or other user defined threshold) of the average read length, then the sequence is marked as a chimeric sequence and sorted into a 'bad' directory. Notably, since reads mapping generally declines towards the ends of the sequence, which causes a natural increase in base disequilibrium, sequence edges are excluded from read disequilibrium checks. Additionally, as the disequilibrium check is only reliable with a sufficient number of mapped reads, seqeunces with a maxium coverage of < 30 are sorted into a 'bad_low' directory to indicate their failure of the disequilibrium filter could be due to low read mapping.

In addition to sorting sequences, ChimeraKiller also prints coverage plots for each sequence making it to the second filtering step so that users can review and modify sequence assignments as they see fit. Coverage plots also identify sequence regions with an identical match to other sequences which can help in the identification of chimeric sequnces. 
