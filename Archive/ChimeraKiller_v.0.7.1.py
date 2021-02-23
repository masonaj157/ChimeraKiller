#!/usr/bin/env python

# ChimeraSorter was made as a followup to ChimeraChecker. It uses several filters to identify and sort filters into good and bad. If you have questions about this script, send a postcard with your information and question addressed to My Ass, You're Screwed Blvd., Imalone, WI 54848

# Additional software necessary to run this:
# (1) bwa 
# (2) gatk
# (3) samtools
# (4) bedtools
# (5) picard
# (6) bam-read from transrate
import time
start = time.time()

import argparse
import sys, os, shutil
import subprocess
import csv
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
import random
from scipy import stats
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
import pysam
from scipy.signal import argrelextrema

# Set up and extract the command line options.
parser = argparse.ArgumentParser(description='Automated checking of coding seqeuences for chimeric transcripts')
parser.add_argument("-i","--input",
                    type=argparse.FileType('r'),
                    help="Fasta file of contigs to check. Use only CODING SEQUENCES.")
parser.add_argument("-r","--reads",
                    type=argparse.FileType('r'),
                    help="Fastq file of UNPAIRED reads.")
parser.add_argument("-b","--bwa",
                    nargs='?',
                    type=str,
                    default="/usr/local/bin/bwa",
                    help="Path to bwa executable. Default assumes it is in your PATH.")
parser.add_argument("-s","--samtools",
                    nargs='?',
                    type=str,
                    default="samtools",
                    help="Path to samtools executable. Default assumes it is in your PATH.")
parser.add_argument("-bt","--bedtools",
                    nargs='?',
                    type=str,
                    default="bedtools",
                    help="Path to bedtools executable. Default assumes it is in your PATH.")
parser.add_argument("-m","--mismatches",
                    nargs='?',
                    type=int,
                    default=0,
                    help="Number of allowable mismatches to keep in alignment. The default is 0.")
parser.add_argument("-ts","--tooShort",
                    nargs='?',
                    type=int,
                    default=150,
                    help="Minimum length of clipped reads. Default is 150.")
parser.add_argument("-mq","--mapQuality",
                    nargs='?',
                    type=int,
                    default=0,
                    help="Minimum mapping quality. Default is 0. Note that reads with multiple mappings get assigned a quality of 0 by bwa.")
parser.add_argument("-min","--minRead",
                    nargs='?',
                    type=int,
                    default=150,
                    help="Minimum read length. Default is 150.")
parser.add_argument("-max","--maxRead",
                    nargs='?',
                    type=int,
                    default=1000,
                    help="Maximum read length. Default is 1000.")
parser.add_argument("-p","--picard",
                    nargs='?',
                    type=str,
                    default="java -jar /path/picard.jar ",
                    help="Picard command. Something like: java -jar /PATH/TO/PICARD/picard.jar")
parser.add_argument("-g","--gatk",
                    nargs='?',
                    type=str,
                    default="java -jar /path/GenomeAnalysisTK.jar ",
                    help="GATK command. Some like: java -jar /PATH/TO/GATK/GenomeAnalysisTK.jar")
parser.add_argument("-d","--readDifference",
					nargs='?',
					type=float,
					default=0.75,
					help="Value between 0 and 1 for the threshold for percentage difference between average lengths of reads on either side of a given site based on the average read size.  Transcripts with a site with a descrepancy between read lengths above this percent are removed based on the assumption that reads should be approximately evenly distributed across sites (example: for a given site, if the average number of bases per read to the left is 10 and the average number of bases per read to the right is 200, the difference between them is 190, which flags the transcript for removal. Default is 0.75")
parser.add_argument("-mer","--merSize",
                    nargs='?',
                    type=int,
                    default=100,
                    help="Mer size for sequence comparisons. Default is 100.")
args = parser.parse_args()



# Check the input fasta file and output the number of detected sequences.
totalContigs = 0
for seq in SeqIO.parse(args.input,format="fasta") :
    totalContigs += 1
print("Total number of input contigs = " + str(totalContigs))
print("Maximum allowed mismatches per read retained in the alignment = " + str(args.mismatches))

print("*"*100)


# Function to generate mers for a single sequence (character string)
def MerSplit(merSize, seq) :
    mers = []
    if len(seq) < merSize : 
        print("Sequence length (" + str(len(seq)) + ") is less than mer size (" + str(merSize) + ").")
        sys.exit()
    for i in range(0,len(seq)-merSize+1) :
        mers.append(seq[i:i+merSize])
    return mers


# Function to identify exactly matching regions between two sequences
# The return value will be a list of three-element lists, where the first element is a list of the starting and ending positions
# in seqence 1, the second is a list for sequence 2, and the third is the actual sequence match.
def ExactMatches(merSize,seq1,seq2) :
    if len(seq1) < merSize or len(seq2) < merSize :
        print("Sequence length is less than mer size (" + str(merSize) + ").")
        sys.exit()
    matches = []
    seq1mers = MerSplit(merSize,seq1)
    seq2mers = MerSplit(merSize,seq2)
    i = 0
    while i < len(seq1mers) :
        if seq1mers[i] in seq2mers : 
            match1 = [i,i+merSize]                # get seq1 coordinates for initial match
            seq = seq1mers[i]                     # get initial match sequence
            j = seq2mers.index(seq1mers[i])       # get seq2 coordinates for initial match
            match2 = [j,j+merSize]
            i += 1
            j += 1
            while i < len(seq1mers) and j < len(seq2mers) and seq1mers[i] == seq2mers[j] :    # extend match while 2 sequences agree
                match1[1] += 1
                match2[1] += 1
                seq += seq1mers[i][-1]
                i += 1
                j += 1
            match1[1] += -1                       # to correct for overstepping by one 
            match2[1] += -1
            matches.append([match1,match2,seq])
        else :
            i += 1
            continue
    return matches


# Run the alignment
grepNM = "grep -E 'NM:i:[0-" + str(args.mismatches) + "][[:space:]]|^@'"
name = args.input.name.split(".")[0] 
print("Generating bwa alignment: " + name + ".bam")
# Build the bwa index
command = args.bwa + " index " + args.input.name
subprocess.call(command,shell=True)
# Generate the initial sam alignment
command = args.bwa + " mem -M -t 4 -R \'@RG\\tID:" + args.input.name + "\\tSM:" + args.reads.name + "' " + args.input.name + " " + args.reads.name + " | " + grepNM  + " > tmp1.sam"
subprocess.call(command,shell=True)
# Create a sorted bam file
command = args.picard + " SortSam INPUT=tmp1.sam OUTPUT=tmp2.bam SORT_ORDER=coordinate" 
subprocess.call(command,shell=True)
command = args.picard + " BuildBamIndex INPUT=tmp2.bam"
subprocess.call(command,shell=True)
# Remove overclipped reads
command = args.picard + " CreateSequenceDictionary REFERENCE=" + args.input.name + " OUTPUT=" + args.input.name.split(".")[0] + ".dict"
subprocess.call(command,shell=True)
command = args.samtools + " faidx " + args.input.name
subprocess.call(command,shell=True)
command = args.gatk + " -T PrintReads -R " + args.input.name + " -I tmp2.bam  -rf OverclippedRead --filter_is_too_short_value " + str(args.tooShort) + " --do_not_require_softclips_both_ends -rf MappingQuality -mmq " + str(args.mapQuality) + " -rf ReadLength -minRead " + str(args.minRead) + " -maxRead " + str(args.maxRead) + " -o " + name + ".bam"
subprocess.call(command,shell=True)
# Calculate coverage
command = args.bedtools + " genomecov -d -ibam " + name + ".bam > coverage.txt"
subprocess.call(command,shell=True)
subprocess.call("rm tmp*[bs]a[im]",shell=True)

# Read in the coverage data
print("*"*100)
print("Importing coverage information.")
X = pd.read_csv("coverage.txt",sep='\t',names=['transcript','pos','cov'])
print("Identified coverage data for " + str(len(set(X['transcript']))) + " transcripts.")
T = list(set(X['transcript']))

print("*"*100)
subprocess.call("mkdir plots", shell=True)
subprocess.call("mkdir plots/good", shell=True)
subprocess.call("mkdir plots/bad", shell=True)
subprocess.call("mkdir plots/bad_low", shell=True)
zeroBad = []
# Check for any sites with no coverage make list of these files
for i in T :
    zeros = list(X[X['transcript'] == i]['cov']).count(0)
    if zeros :
        zeroBad.append(i)
        
# Remove the contigs with zero coverage sites from the master list
T = [x for x in T if x not in zeroBad]        
 

# Create a dictionary of all of the contigs.
seqFile = open(args.input.name,"r")
S = {}
for seq in SeqIO.parse(seqFile,"fasta") :
    S[seq.description] = str(seq.seq)
seqFile.close()


## Read in samfile
samname = name + ".bam"    
samfile = pysam.AlignmentFile(samname, "rb")        
read_lengths = []
references = []
for read in samfile:
	references.append(read.reference_name)        
	read_lengths.append(read.rlen)
	

read_mean = np.mean(read_lengths)
Diff_cutoff = read_mean * args.readDifference
print("Average read size is " + str(read_mean))
print("Looking for sites with difference in average read distribution greater than " + str(Diff_cutoff))
print("*"*100)

	
print("Starting transcript filtering")
print("*"*100)

##This stuff is the working (as in in development) code.
good_transcripts = []
bad_transcripts = [] 
for transcript in T:
	print(transcript)
	x = list(X[X['transcript'] == transcript]['pos'])
	x = np.asarray(x)
	ydata = list(X[X['transcript'] == transcript]['cov'])
	ydata = np.asarray(ydata)
	matches = {}
	for j in T :
		match = ExactMatches(100,S[transcript],S[j])
		#match = ExactMatches(args.merSize,S[transcript],S[j])
		#print(match)
		if match and S[transcript] != S[j] :
			matches[j] = match
                 
	if len(matches.keys()) == 0 :
		plt.figure(figsize=(12,4))
		plt.title(transcript)
		plt.ylabel('Coverage')
		plt.xlabel('CDS position')
		plt.bar(x,ydata,color="blue",ec="blue")
		plt.xlim(x[0],x[-1])
	if not len(matches.keys()) :
		plt.title(transcript)
		plt.ylabel('Coverage')
		plt.xlabel('CDS position')
		plt.bar(x,ydata,color="blue",ec="blue")
		plt.xlim(x[0],x[-1])
	else :
		plt.figure(figsize=(12,4+4))
		plt.subplot(2,1,1)
		plt.title(transcript)
		plt.ylabel('Coverage')
		plt.xlabel('CDS position')
		plt.bar(x,ydata,color="blue",ec="blue")
		plt.xlim(x[0],x[-1])
		plt.subplot(2,1,2)
		plt.xlim(x[0],x[-1])
		plt.ylim(0,len(matches.keys())+1)
		pos = 1
		for m in matches.keys() :
			plt.text(0.1*x[-1],pos+0.25,m)
			for n in range(0, len(matches[m])) :
				a = [matches[m][n][0][0], matches[m][n][0][1]]
				plt.hlines(pos,a[0],a[1],"r",linewidth=3)
				plt.vlines(a[0],pos+0.1,pos-0.1)
				plt.vlines(a[1],pos+0.1,pos-0.1)
			pos += 1
	lengths = []
	left = []
	right = []
	position = []
	read = []
	mean_left = []
	mean_right = []
	#transcript = 'Bnigr-CLP1856_SVMPIII-1_extContig424_-_CDS'
	for i in range(0,samfile.get_reference_length(transcript)):               
		iter = samfile.fetch(transcript, i, i+1)
		samples = []
		tmp_left = []
		tmp_right = []
		for x in iter:
			samples.append(x)
		## Downsample to avoid huge calculations and maybe biases?
		if len(samples) > 100:
			new_samples = random.sample(samples, 100)
		else:
			new_samples = samples
		for x in samples:
			lengths.append(int(x.reference_end-x.reference_start))
			left.append(int(i - x.reference_start))
			tmp_left.append(int(i - x.reference_start))
			right.append(int(x.reference_end-i))
			tmp_right.append(int(x.reference_end-i))
			read.append(x.qname)
		mean_left.append(np.mean(tmp_left))
		mean_right.append(np.mean(tmp_right))          
		b = [i] * len(new_samples)
		position.extend(b)
	single_positions = []
	single_positions = range(0,samfile.get_reference_length(transcript))
	Diff = abs(np.asarray(mean_left)-np.asarray(mean_right))
	df = pd.DataFrame(Diff, columns=['Diff'])
	n = 50
	df['min'] = df.iloc[argrelextrema(df.Diff.values, np.less_equal, order=n)[0]]['Diff']
	df['max'] = df.iloc[argrelextrema(df.Diff.values, np.greater_equal, order=n)[0]]['Diff']
	min = df.iloc[argrelextrema(df.Diff.values, np.less_equal, order=n)[0]]['Diff']
	max = df.iloc[argrelextrema(df.Diff.values, np.greater_equal, order=n)[0]]['Diff']
	subsetDiff = df.Diff[min.index[0]:min.index[len(min.index)-1]]
	subsetDiff = pd.DataFrame({'Diff':subsetDiff.values})
	if	np.min(ydata) < np.median(ydata)*0.0025:
		print('bad')
		plt.savefig("plots/bad/" + transcript +".pdf")
		plt.close()
		bad_transcripts.append(transcript)
	elif subsetDiff['Diff'].max() > Diff_cutoff and np.max(ydata) > 30:
		print('bad')
		plt.savefig("plots/bad/" + transcript +".pdf")
		plt.close()
		bad_transcripts.append(transcript)
	elif subsetDiff['Diff'].max() > Diff_cutoff and np.max(ydata) < 30:
		print('low')
		plt.savefig("plots/bad_low/" + transcript +".pdf")
		plt.close()
		bad_transcripts.append(transcript)
	else:
		print('good')
		plt.savefig("plots/good/" + transcript +".pdf")
		plt.close()
		good_transcripts.append(transcript)
		
		
print(str(len(good_transcripts)) + " sequences being written to good.fasta")
print(str(len(bad_transcripts)) + " sequences being written to bad.fasta")

sequences = list(SeqIO.parse(args.input.name,"fasta"))

good_seqs=[]
for seq in sequences:
	for good in good_transcripts:
		if seq.id == good:
			good_seqs.append(seq)


bad_seqs=[]
for seq in sequences:
	for bad in bad_transcripts:
		if seq.id == bad:
			bad_seqs.append(seq)
			

print("*"*100)
subprocess.call("mkdir fastas", shell=True)
subprocess.call("mkdir fastas/good", shell=True)
subprocess.call("mkdir fastas/bad", shell=True)

for seq in good_seqs:
	record =[]
	record.append(seq)
	good_fasta = "fastas/good/" + seq.id + ".fasta"
	handle=open(good_fasta, "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(record)
	handle.close()
	
	
for seq in bad_seqs:
	record =[]
	record.append(seq)
	bad_fasta = "fastas/bad/" + seq.id + ".fasta"
	handle=open(bad_fasta, "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(record)
	handle.close()
 
end = time.time()
print(end - start)
                