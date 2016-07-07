#!/bin/bash

breseqFolder=$1
refpath=$2
mutationList=$3

# Get list of bam files to use and get strain names
cd $breseqFolder/postProcess/08pileup
find ../../*/data/reference.bam > bamlist
cut -d'/' -f3 < bamlist > strainlist

# Create list of positions for samtools
# Bed format is 0-based, half-open, unlike sam format. Thus we subtract by 1 here.
chr=NC_007795
awk -v chr="$chr" -F"," '{ print chr " " $1-1 " " $1}' mutationlist_20160528.txt > positions

# Call mpileup for each mutation in the list
samtools mpileup --max-depth 10000 -f $refpath --no-BAQ --positions positions --bam-list bamlist > mutation_pileup

# C++ scipt to parse mpileup output
# This outputs two files:
# forcecalled_counts.txt
# forcecalled_frequencies.txt
# For each file, a row corresponds to a single mutation in the mutationlist file
# Each column corresponds to strainlist
./parsePileup strainlist mutation_pileup $mutationList