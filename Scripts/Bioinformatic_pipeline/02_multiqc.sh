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

## Create multiqc directory, you need to create it before running multiqc
mkdir -p "/Users/martaibanez/Desktop/PhD/NXT06/FASTQS/HEV/results/multiqc"

## put path to multiqc output directory
MULTIQC_OUT="/Volumes/MARTA_06/PhD/RUNS/NXT15/results/multiqc"

## change to fastqc directory to run multiqc
cd "/Volumes/MARTA_06/PhD/RUNS/NXT15/results/fastqc"

## run multiqc for all fastqs in the folder: -o states the output directory, which is the one created before
multiqc --module fastqc .  -o $MULTIQC_OUT
