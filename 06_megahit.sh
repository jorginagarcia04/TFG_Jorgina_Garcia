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
#
## --- IMPORTANT !!!! ---
## This script needs to be run with UNZIPPED FASTQ FILES (.fastq)
## megahit does not accept .fastq.gz files as input

## set path to fastq files  after removal of human genome (bowtie2)
FILTERED_FASTQ_PATH="/Users/martaibanez/Desktop/Metagenomics/Gambusia_Nov17/kraken2_filtered_Nov18"


## create output directory for assemblies
##mkdir -p "/Users/martaibanez/Desktop/PhD/NXT06/FASTQS/HEV2/results/assembly"
## set path to output directory for assemblies
ASSEMBLY_PATH="/Users/martaibanez/Desktop/Metagenomics/Gambusia_Nov17/megahit_filtered"

## get unique file names without R1/R2 designation
FILE_NAMES=$(find $FILTERED_FASTQ_PATH -type f -name "*.fastq" | rev |cut -c 9- | rev | sort -u)
## print file names to check
echo $FILE_NAMES

## start loop to run megahit for each pair of fastq files
for sample in $FILE_NAMES
do
    ## set input names
    in1="${sample}R1.fastq"
    in2="${sample}R2.fastq"
    ## set output path
    out="$ASSEMBLY_PATH/$(basename $sample)out"
    ## print output
    echo $in1 $in2 $out
    ## run megahit
    megahit -m 0.6 -t 5 -1 $in1 -2 $in2 -t 6 -o  $out
done




