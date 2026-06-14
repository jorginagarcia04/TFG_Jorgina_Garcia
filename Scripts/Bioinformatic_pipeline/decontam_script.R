## --- DECONTAMINATION OF SAMPLES --- ##
## USING A PACKAGE NAMED DECONTAM ##
## R version 4.5.2 (2025-10-31)
## Copyright: Marta Ibañez Lligoña (marta.ibanez@vhir.org)
## Adapted by: Jorgina Garcia Larroy


## --- Install needed packages (if not already installed) --- ##

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")
install.packages("readxl")

## --- Load packages --- ##
library(decontam)
library(readxl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)

`%nin%` = Negate(`%in%`)

## --- INPUT SAMPLES METADATA --- ##
## specify the codes and samples found
samples.metadata <- read_xlsx("C:/Users/jorgi/OneDrive - Universitat Pompeu Fabra/BIO HUMANA/4t BIO HUMANA/TFG/Analysis/Samples_class.xlsx")

## --- Set working directory --- ##
setwd("C:/Users/jorgi/OneDrive - Universitat Pompeu Fabra/BIO HUMANA/4t BIO HUMANA/TFG/Analysis/Kraken")

## Now we will define where the taxonomy tables will be imported from
## Specify directory where the tables are found (READ LEVEL, SHORT SEQS) and store it in a variable
kraken.dir <- "C:/Users/jorgi/OneDrive - Universitat Pompeu Fabra/BIO HUMANA/4t BIO HUMANA/TFG/Analysis/Kraken/kraken_reads"

## list fastp tables (all the files in the directory finishing in _kreport.txt)
kreport_files <- list.files(path = kraken.dir,pattern = "_kreport.txt") ## XX files

## --- LOAD TABLES --- ##

## Function to load all the tables into R. You have to call this function and specify the rank that you're working with. 
load_kraken_reports <- function(kreport_files, kraken_dir, rank = "G") {
  # Move to directory where the files are
  setwd(kraken.dir)
  
  # Read all files as a table and save into a named list
  data_list <- lapply(kreport_files, function(file) {
    read.delim2(file, header = FALSE)
  })
  # Name the list elements as the file name
  names(data_list) <- kreport_files
  
  # Process each table
  data_list <- lapply(data_list, function(df) {
    
    # Add column names
    colnames(df) <- c("Percentage", "Nreads", "UNreads", "Rank", "NCBI_ID", "Name")
    
    # Clean white space in taxon name
    df$Name <- trimws(df$Name)
    
    # Filter by rank (e.g., G = genus, F = family)
    df <- df[df$Rank == rank, ]
    
    # Remove unnecessary columns
    df <- df[, c("Nreads", "Rank", "NCBI_ID", "Name")]
    
    # Rename last column depending on rank (e.g. G to Genus)
    colnames(df)[4] <- switch(rank,
                              "G" = "GENUS",
                              "F" = "FAMILY",
                              "S" = "SPECIES",
                              "TAXON")
    
    return(df)
  })
# Returns the data list of the processed tables  
  return(data_list)
}

## Use the function to load the tables and save the results
kraken_results <- load_kraken_reports(kreport_files, kraken.dir, rank = "G")
## This needs to be done to put the table into the environment so you can work with these. Theoretically you could also work with the kraken_results objects (it works like a list)
list2env(kraken_results, envir = .GlobalEnv)

## --- STARTING DECONTAMINATION ANALYSIS --- ##

## list all files that are found in the environment (previously loaded)
list_files.ALL <- (ls(pattern = "_kreport.txt"))

## get the data frames into a list of objects that you can work with
list_files.ALL <- mget(list_files.ALL)

## check number of genera for each control, just for information
lapply(list_files.ALL, nrow)

## Join all rows of the different tables to then make a matrix of counts.
## The result is a table with genera as rows and columns counts per sample
full.kraken_results <- list_files.ALL %>% reduce(full_join, by= c('GENUS', 'NCBI_ID'))

## Remove rank columns
full.kraken_results <- full.kraken_results %>%
  select(-starts_with("Rank"))

## change first 3 column's names
colnames(full.kraken_results)[1:3] <- c("C_EXT", "ID_NCBI", "NAME")

## remove columns that start with NCBI_ID that appear after the joining
full.kraken_results <- full.kraken_results %>%
  select(-starts_with("NCBI_ID."))

## change column names to the actual names of the samples
colnames(full.kraken_results)[4:ncol(full.kraken_results)] <- names(list_files.ALL[2:length(list_files.ALL)])

## relocate the first column
full.kraken_results <- full.kraken_results %>% relocate(C_EXT, .after = NAME)

# Replace NAs with 0 for the whole dataframe
full.kraken_results <- full.kraken_results %>%
  replace(is.na(.), 0)

# Create new matrix where the taxons are row names so these don't interfere with the analysis
full.kraken_mat <- as.matrix(full.kraken_results[, -c(1,2)])
rownames(full.kraken_mat) <- full.kraken_results[,2] ## this can be done with the name of the genus/family or the ncbi id, which is more correct

## Change names of columns to the NCB IDs
mat_ids <- sub("_S.*", "", colnames(full.kraken_mat)) ## CHANGE ACCORDING TO LEVEL OF ANALYSIS
samples.metadata <- samples.metadata[match(mat_ids, samples.metadata$Sample_ID),]
colnames(full.kraken_mat) <- mat_ids

## Now we can save the matrix in csv format in case we need it at some point 
write.csv(full.kraken_mat,"read_counts.csv")

## --- Now we can actually start using decontam --- ##

## Eliminates columns where lecture total is 0 
counts <- full.kraken_mat[,colSums(full.kraken_mat) > 0]

## samples metadata reformatting
samples.metadata <- as.data.frame(samples.metadata)

## set rownames to sample metadata
rownames(samples.metadata) <- samples.metadata$Sample_ID

## only keep those samples that passed the filter
samples.metadata <- samples.metadata[colnames(counts),]

## You can run this to check everything makes sense and all samples are where these should be
all(colnames(counts) %in% rownames(samples.metadata))  # Should be TRUE
all(rownames(samples.metadata) %in% colnames(counts))  # Also should be TRUE

## 1 Let's run first step, use the sample_data function to change to the required format
samp <- sample_data(samples.metadata)
rownames(samp) <- samp$Sample_ID

# 2. Create phyloseq object (decontam prefers this format), taxons are the rows
ps <- phyloseq(otu_table(counts, taxa_are_rows=TRUE), samp)

# 3. Identify which are controls
sample_data(ps)$is.neg <- sample_data(ps)$Type == "Control"

# 4. Run prevalence-based contamination detection, comparing prevalence between controls and samples
contam_results <- isContaminant(ps, method = "prevalence", neg = "is.neg")

# 5. Extract results
contaminants <- rownames(subset(contam_results, contaminant == TRUE))
## Final table of contaminants. This is what should be removed from the final tables. 
contaminants.table <- subset(contam_results, contaminant == TRUE)

## save table of contaminants
write.csv(contaminants.table, "decontam_contaminants.csv")


## --- REMOVAL OF THE IDENTIFIED CONTAMINANTS FROM EACH SAMPLE --- ##

## --- Removal of contaminants PP7_ENDOMETRIUM

#Load table and filters the genera in the contaminant list
PP7_E <- PP7_E_S2_L001_fastp__kreport.txt[which(PP7_E_S2_L001_fastp__kreport.txt$GENUS %nin% contaminants),]
#Eliminate genera with less than 5 reads (noise elimination)
PP7_E <- PP7_E[PP7_E$Nreads >=5,]

#CPM calculation: reads for each genus/total reads * 10^6
PP7_E <- PP7_E %>%
  mutate(CPM = (Nreads/33722)*1000000) 
#Eliminate the reads corresponding to human (Homo genus)
PP7_E <- PP7_E[PP7_E$GENUS != "Homo", ]
#Change rownames into the genus names and delete the genus column
rownames(PP7_E) <- PP7_E$GENUS
PP7_E$GENUS <- NULL

#Save filtered .csv
write.csv(PP7_E,"PP7_E_net2.csv")

## --- Removal of contaminants PP7_PM

#Load table and filters the genera in the contaminant list
PP7_PM <- PP7_PM_S3_L001_fastp__kreport.txt[which(PP7_PM_S3_L001_fastp__kreport.txt$GENUS %nin% contaminants),]
#Eliminate genera with less than 5 reads (noise elimination)
PP7_PM <- PP7_PM[PP7_PM$Nreads >=5,]

#CPM calculation: reads for each genus/total reads * 10^6
PP7_PM <- PP7_PM %>%
  mutate(CPM = (Nreads/38835)*1000000)
#Eliminate the reads corresponding to human (Homo genus)
PP7_PM <- PP7_PM[PP7_PM$GENUS != "Homo", ]
#Change rownames into the genus names and delete the genus column
rownames(PP7_PM) <- PP7_PM$GENUS
PP7_PM$GENUS <- NULL

#Save filtered .csv
write.csv(PP7_PM,"PP7_PM_net2.csv")

## --- Removal of contaminants PP7_PT

#Load table and filters the genera in the contaminant list
PP7_PT <- PP7_PT_S1_L001_fastp__kreport.txt[which(PP7_PT_S1_L001_fastp__kreport.txt$GENUS %nin% contaminants),]
#Eliminate genera with less than 5 reads (noise elimination)
PP7_PT <- PP7_PT[PP7_PT$Nreads >=5,]

#CPM calculation: reads for each genus/total reads * 10^6
PP7_PT <- PP7_PT %>%
  mutate(CPM = (Nreads/51556)*1000000)
#Eliminate the reads corresponding to human (Homo genus)
PP7_PT <- PP7_PT[PP7_PT$GENUS != "Homo", ]
#Change rownames into the genus names and delete the genus column
rownames(PP7_PT) <- PP7_PT$GENUS
PP7_PT$GENUS <- NULL

#Save filtered .csv
write.csv(PP7_PT,"PP7_PT_net2.csv")

## --- Removal of contaminants PP7_V

#Load table and filters the genera in the contaminant list
PP7_V <- PP7_V_S4_L001_fastp__kreport.txt[which(PP7_V_S4_L001_fastp__kreport.txt$GENUS %nin% contaminants),]
#Eliminate genera with less than 5 reads (noise elimination)
PP7_V <- PP7_V[PP7_V$Nreads >=5,]

#CPM calculation: reads for each genus/total reads * 10^6
PP7_V <- PP7_V %>%
  mutate(CPM = (Nreads/57243)*1000000)
#Eliminate the reads corresponding to human (Homo genus)
PP7_V <- PP7_V[PP7_V$GENUS != "Homo", ]
#Change rownames into the genus names and delete the genus column
rownames(PP7_V) <- PP7_V$GENUS
PP7_V$GENUS <- NULL

#Save filtered .csv
write.csv(PP7_V,"PP7_V_net2.csv")

## --- Removal of contaminants PP8_C

#Load table and filters the genera in the contaminant list
PP8_C <- PP8_C_S7_L001_fastp__kreport.txt[which(PP8_C_S7_L001_fastp__kreport.txt$GENUS %nin% contaminants),]
#Eliminate genera with less than 5 reads (noise elimination)
PP8_C <- PP8_C[PP8_C$Nreads >=5,]

#CPM calculation: reads for each genus/total reads * 10^6
PP8_C <- PP8_C %>%
  mutate(CPM = (Nreads/701644)*1000000)
#Eliminate the reads corresponding to human (Homo genus)
PP8_C <- PP8_C[PP8_C$GENUS != "Homo", ]
#Change rownames into the genus names and delete the genus column
rownames(PP8_C) <- PP8_C$GENUS
PP8_C$GENUS <- NULL

#Save filtered .csv
write.csv(PP8_C,"PP8_C_net2.csv")

## --- Removal of contaminants PP8_F

#Load table and filters the genera in the contaminant list
PP8_F <- PP8_F_S9_L001_fastp__kreport.txt[which(PP8_F_S9_L001_fastp__kreport.txt$GENUS %nin% contaminants),]
#Eliminate genera with less than 5 reads (noise elimination)
PP8_F <- PP8_F[PP8_F$Nreads >=5,]

#CPM calculation: reads for each genus/total reads * 10^6
PP8_F <- PP8_F %>%
  mutate(CPM = (Nreads/7235)*1000000)
#Eliminate the reads corresponding to human (Homo genus)
PP8_F <- PP8_F[PP8_F$GENUS != "Homo", ]
#Change rownames into the genus names and delete the genus column
rownames(PP8_F) <- PP8_F$GENUS
PP8_F$GENUS <- NULL

#Save filtered .csv
write.csv(PP8_F,"PP8_F_net2.csv")

## --- Removal of contaminants PP8_PM

#Load table and filters the genera in the contaminant list
PP8_PM <- PP8_PM_S6_L001_fastp__kreport.txt[which(PP8_PM_S6_L001_fastp__kreport.txt$GENUS %nin% contaminants),]
#Eliminate genera with less than 5 reads (noise elimination)
PP8_PM <- PP8_PM[PP8_PM$Nreads >=5,]

#CPM calculation: reads for each genus/total reads * 10^6
PP8_PM <- PP8_PM %>%
  mutate(CPM = (Nreads/95873)*1000000)
#Eliminate the reads corresponding to human (Homo genus)
PP8_PM <- PP8_PM[PP8_PM$GENUS != "Homo", ]
#Change rownames into the genus names and delete the genus column
rownames(PP8_PM) <- PP8_PM$GENUS
PP8_PM$GENUS <- NULL

#Save filtered .csv
write.csv(PP8_PM,"PP8_PM_net2.csv")

## --- Removal of contaminants PP8_PT

#Load table and filters the genera in the contaminant list
PP8_PT <- PP8_PT_S8_L001_fastp__kreport.txt[which(PP8_PT_S8_L001_fastp__kreport.txt$GENUS %nin% contaminants),]
#Eliminate genera with less than 5 reads (noise elimination)
PP8_PT <- PP8_PT[PP8_PT$Nreads >=5,]

#CPM calculation: reads for each genus/total reads * 10^6
PP8_PT <- PP8_PT %>%
  mutate(CPM = (Nreads/4333795)*1000000)
#Eliminate the reads corresponding to human (Homo genus)
PP8_PT <- PP8_PT[PP8_PT$GENUS != "Homo", ]
#Change rownames into the genus names and delete the genus column
rownames(PP8_PT) <- PP8_PT$GENUS
PP8_PT$GENUS <- NULL

#Save filtered .csv
write.csv(PP8_PT,"PP8_PT_net2.csv")

