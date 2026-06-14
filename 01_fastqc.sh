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

## Create directory to save the report files
## mkdir -p "/Users/martaibanez/Desktop/PhD/NXT06/FASTQS/fastqc"
## set path to new directory as cleanfastq
cleanfastq="/home/sumeyye/Desktop/P1_26/results/fastqc"

## run fastqc through all documents with just one command, you need to add the fastqs path where your files are located, *.fastq.gz is selecting all files at once.
fastqc /home/sumeyye/Desktop/P1_26/data/*.fastq.gz -t 8 -o  $cleanfastq