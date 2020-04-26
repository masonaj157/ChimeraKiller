# ChimeraKiller
In-development script for identifying chimeric transcripts in a de novo assembled transcriptome


Suggestion for creating conda environment:


conda create -n chimerakiller_env python=3.6 biopython bwa gatk samtools bedtools picard transrate-tools joblib pandas matplotlib scipy pysam gatk4

conda activate chimerakiller_env

pip install pathos
