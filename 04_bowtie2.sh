#!/bin/bash
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


## VARIABLES/PATHS NEEDED FOR THIS SCRIPT
## CLEAN_PATH-> path to clean FASTQ files
## REF_GENOME -> path where the indexed files for the reference genome are located
## MAPPED_PATH --> where the mapped bam files will be saved. 
## UNMAPPED_FASTQ --> where the unmapped bam files will be saved


## Indicate folder with clean paired-end fastq reads
CLEAN_PATH="/Volumes/MARTA_06/PhD/Optimitzacio/VHIR131/results/cont_out"

## Create directory for log files for bowtie2
mkdir -p "/Users/martaibanez/Desktop/PhD/NXT07/results/logs/bowtie2"
## establish path for log files
LOG_PATH="/Volumes/MARTA_06/PhD/Optimitzacio/VHIR131/results/logs/bowtie2"

## Retrieve file names
FILE_NAMES=$(find $CLEAN_PATH -type f -name "*.fastq.gz" | rev |cut -c 13- | rev | sort -u)
## print file names into the screen
echo These are the samples $FILE_NAMES

## This is the path for the indexed human genome, this needs to be the reference genome, whichever is needed for the fastq
REF_GENOME="/Volumes/MARTA_06/metagenomics_data/human/Homo_sapiens"

## Create folder for BAM files
mkdir -p "/Users/martaibanez/Desktop/PhD/NXT06/FASTQS/HEV2/results/mapped"
## establish directory
MAPPED_PATH="/Volumes/MARTA_06/PhD/Optimitzacio/VHIR131/results/mapped"

## Create folder for BAM files
mkdir -p "/Users/martaibanez/Desktop/PhD/NXT06/FASTQS/HEV2/results/mapped_sorted"
## establish directory
SORTED_PATH="/Volumes/MARTA_06/PhD/Optimitzacio/VHIR131/results/mapped_sorted"

## Create folder for BAM files with unmapped reads
mkdir -p "/Users/martaibanez/Desktop/PhD/NXT06/FASTQS/HEV2/results/unmapped"
## establish direcotry
UNMAPPED_PATH="/Volumes/MARTA_06/PhD/Optimitzacio/VHIR131/results/unmapped"

## Create folder for new filtered FASTQ files 
mkdir -p "/Users/martaibanez/Desktop/PhD/NXT06/FASTQS/HEV2/results/filtered_fastqs"
## establish directory
UNMAPPED_FASTQ="/Volumes/MARTA_06/PhD/Optimitzacio/VHIR131/results/filtered_fastqs"

## loop through fastq files
for sample in $FILE_NAMES
do
    ## set input file names
    echo Mapping starting for $sample
    in1="${sample}_R1.fastq.gz"
    in2="${sample}_R2.fastq.gz"
    ## Alignment to reference with bowtie2
    bowtie2 --very-sensitive-local --mm --threads 6 --reorder --seed 4 -x $REF_GENOME -1 $in1 -2 $in2 -S "$MAPPED_PATH/$(basename $sample).sam" #> "$LOG_PATH/$(basename $sample)_alignment.log"
    ## Convert SAM file to BAM file
    samtools view -@ 2 -bS "$MAPPED_PATH/$(basename $sample).sam" > "$MAPPED_PATH/$(basename $sample).bam"
    ## | samtools view -o "$MAPPED_PATH/$(basename $sample).bam"
    echo Mapping to reference genome for $in1 and $in2 done
done

## Move to mapped folder to sort bam files
echo Moving to folder $MAPPED_PATH
MAPPED_NAMES=$(find $MAPPED_PATH -type f -name "*.bam" | rev |cut -c 5- | rev | sort -u)
echo $MAPPED_NAMES

## Loop to sort bam files and extract unmapped reads into new fastq files
for file in $MAPPED_NAMES
do
    ## SORT BAM FILE
    echo The sample being sorted is "$MAPPED_PATH/$(basename $file).bam"
    samtools sort -n "$MAPPED_PATH/$(basename $file).bam" -o "$SORTED_PATH/$(basename $file)_sorted.bam"
    ## -f pot ser o 5 o 8
    ## EXTRACT UNMAPPED READS FROM SORTED BAM FILE
    samtools view -@ 4 -b -f 12 -F 256   "$SORTED_PATH/$(basename $file)_sorted.bam" > "$UNMAPPED_PATH/$(basename $file)_UNMAP.bam"
    out1="$UNMAPPED_FASTQ/$(basename $file)_R1.fastq.gz"
    out2="$UNMAPPED_FASTQ/$(basename $file)_R2.fastq.gz"
    ## Transform to fastq files 
    samtools fastq  -@ 6  -1 $out1 -2 $out2 -n "$UNMAPPED_PATH/$(basename $file)_UNMAP.bam"
    echo BAM file with unmapped reads created with sample $file
done