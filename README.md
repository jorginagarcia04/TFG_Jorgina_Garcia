# Title TFG

### By: Jorgina Garcia Larroy

### Directed by: Dr Josep Quer Sivila

## Repository content
### Readme: 
Information about repository and folders content.

## Scripts ##
### Bioinformatics pipeline ###
This folder contains all scripts used for the analysis of fastQ files from metagenomics samples.

| Script file | Description | 
| ----------------------------- | ----------------------- | 
01_fastqc | Script used to evaluate the quality of raw sequencing reads.   | 
02_multiqc | Script used to evaluate quality across all samples and generate a quality control report. |
03_fastp | Script used to remove low quality, short reads and adapter sequences from fastQ files. | 
04_bowtie2 | Script used to identify human reads by mapping against the human reference genome. |
05_kraken2 | Script used to taxonomically classify the reads. |
06_megahit | Script used to perform de novo assembly of reads into contigs. |
05_kraken2 | Script used to taxonomically classify the resulting contigs. |
decontam_script | Script used to identify potential contaminants across the samples. |

## Interpretation of results ##
This folder contains all scripts used for analysis of results and making of graphs.

| Script file | Description | 
| ----------------------------- | ----------------------- | 
taxonomy_analysis.R | Script used to analyse results from kaiju of main metagenomic analysis. | 
abundance_graphs.R | Script used to make graphs to visually represent the results. | 
pipeline_graphs.R | Scripts used to make graphs to show effectiveness of the bioinformatics pipeline. |
