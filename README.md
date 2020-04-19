# ChimeraKiller
In-development script for identifying chimeric transcripts in a de novo assembled transcriptome


Suggestion for creating conda environment:


conda create -n chimerakiller_env python=3.6 biopython bwa gatk samtools bedtools picard transrate-tools joblib pandas matplotlib scipy pysam

conda activate chimerakiller_env

wget https://console.cloud.google.com/storage/browser/_details/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
gatk3-register path/to/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
