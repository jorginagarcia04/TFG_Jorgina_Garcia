#!/bin/bash
# Copyright (c) 2025 Marta Ibañez Lligoña
# All rights reserved.
#
# This script is part of an internal research pipeline developed at
# Liver Diseases - Viral Hepatitis Lab from Vall d'Hebron Institute of Research.
#
# Unauthorized copying, modification, distribution, or use of this
# script, in whole or in part, is strictly prohibited without
# explicit written permission from the author.
#
# Author: Marta Ibañez Lligoña
# Contact: marta.ibanez@vhir.org

## set path to fastq files  that have had the human reads removed (bowtie2)
FILTERED_FASTQ_PATH="/Volumes/MARTA_06/PhD/RUNS/NXT14/results/filtered_fastqs"

## create output directory for kraken2 results
mkdir -p "/Users/martaibanez/Desktop/PhD/NXT06/FASTQS/HEV2/results/KRAKEN_READS"
## set path to output directory for kraken2 results
KRAKEN_PATH="/Volumes/MARTA_06/PhD/RUNS/NXT14/results/KRAKEN_READS"
 
## Set file names to be analysed, this command must be changed according to the format of the names
FILE_NAMES=$(find $FILTERED_FASTQ_PATH -type f -name "*.fastq.gz" | rev |cut -c 12- | rev | sort -u)
## print names to check
echo $FILE_NAMES

## Set path KRAKEN2 DB, last update 2025, this needs to change according to the location of the DB on your system
KRAKEN2_DB="/Volumes/MARTA_06/metagenomics_data/k2_pluspf_20241228"

## loop through each sample and run kraken2
for sample in $FILE_NAMES
do
    in1="${sample}R1.fastq.gz"
    in2="${sample}R2.fastq.gz"
    ## This is the name for the kraken.txt output file 1 (hard to read)
    out="$KRAKEN_PATH/$(basename $sample)_kraken.txt"
    ## this is the name for the files that will contain the classified reads
    seqsC="$KRAKEN_PATH/$(basename $sample)_classified#.fastq.gz"
    ## this is the name for the files that will contain the unclassified reads
    seqsU="$KRAKEN_PATH/$(basename $sample)_unclassified#.fastq.gz"
    ## this is the name for the kraken report file (easier to read)
    report="$KRAKEN_PATH/$(basename $sample)_kreport.txt"
    echo Running KRAKEN with $in1 and $in2
    kraken2 --db $KRAKEN2_DB --threads 6 --paired $in1 $in2  --report $report --output $out ##--unclassified-out $seqsU ##--classified-out $seqsC  ##--confidence 0.2
    echo KRAKEN2 COMPLETED!
done