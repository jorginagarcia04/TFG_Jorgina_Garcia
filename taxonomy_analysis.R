## --- TAXONOMY ANALYSIS EXAMPLE SCRIPT --- ##
## IDENTIFYING CONTAMINANTS IN METAGENOMICS AT FAMILY/GENUS LEVEL ##
## R version 4.5.2 (2025-10-31)
## Copyright: Marta Ibañez Lligoña (marta.ibanez@vhir.org)
## Adapted by: Jorgina Garcia Larroy

## --- Set working directory --- ##
setwd("C:/Users/jorgi/OneDrive - Universitat Pompeu Fabra/BIO HUMANA/4t BIO HUMANA/TFG/Analysis/Kraken/kraken_reads")

## -- CALCULATING RELATIVE ABUNDANCES of SAMPLES AFTER DECONTAM

#Read csv file and calculate total number of reads.
vagina7 <- read.csv("PP7_V_net2.csv")
sumareads <- sum(vagina7$Nreads)
#Calculate relative abundance
vagina7$abundance <- (vagina7$`Nreads`/sumareads)* 100
#Set genus names as row names
rownames(vagina7) <- vagina7$X
#Remove unnecessary columns
vagina7$X <- NULL
vagina7$Rank <- NULL
#Order the table by decreasing abundance
vagina7 <- vagina7[order(vagina7$abundance, decreasing = TRUE), ]
#Save the table for further analysis
write.csv2(vagina7,"PP7_vagina_abundance.csv")

# Read the csv file for the placental sample of PP7 and calculate relative abundance
placenta7 <- read.csv("PP7_PT_net2.csv")
sumareads <- sum(placenta7$Nreads)
placenta7$abundance <- (placenta7$`Nreads`/sumareads)* 100
rownames(placenta7) <- placenta7$X
placenta7$X <- NULL
placenta7$Rank <- NULL
placenta7 <- placenta7[order(placenta7$abundance, decreasing = TRUE), ]
write.csv2(placenta7,"PP7_placenta_abundance.csv")

# Read the csv file for the plasma sample of PP7 and calculate relative abundance
plasma7 <- read.csv("PP7_PM_net2.csv")
sumareads <- sum(plasma7$Nreads)
plasma7$abundance <- (plasma7$`Nreads`/sumareads)* 100
rownames(plasma7) <- plasma7$X
plasma7$X <- NULL
plasma7$Rank <- NULL
plasma7 <- plasma7[order(plasma7$abundance, decreasing = TRUE), ]
write.csv2(plasma7,"PP7_plasma_abundance.csv")

# Read the csv file for the endometrial sample of PP7 and calculate relative abundance
endometri7 <- read.csv("PP7_E_net2.csv")
sumareads <- sum(endometri7$Nreads)
endometri7$abundance <- (endometri7$`Nreads`/sumareads)* 100
rownames(endometri7) <- endometri7$X
endometri7$X <- NULL
endometri7$Rank <- NULL
endometri7 <- endometri7[order(endometri7$abundance, decreasing = TRUE), ]
write.csv2(endometri7,"PP7_endometri_abundance.csv")

# Read the csv file for the placental sample of PP8 and calculate relative abundance
placenta8 <-read.csv("PP8_PT_net2.csv")
sumareads <- sum(placenta8$Nreads)
placenta8$abundance <- (placenta8$`Nreads`/sumareads)* 100
rownames(placenta8) <- placenta8$X
placenta8$X <- NULL
placenta8$Rank <- NULL
placenta8 <- placenta8[order(placenta8$abundance, decreasing = TRUE), ]
write.csv2(placenta8,"PP8_placenta_abundance.csv")

# Read the csv file for the plasma sample of PP8 and calculate relative abundance
plasma8 <- read.csv("PP8_PM_net2.csv")
sumareads <- sum(plasma8$Nreads)
plasma8$abundance <- (plasma8$`Nreads`/sumareads)* 100
rownames(plasma8) <- plasma8$X
plasma8$X <- NULL
plasma8$Rank <- NULL
plasma8 <- plasma8[order(plasma8$abundance, decreasing = TRUE), ]
write.csv2(plasma8,"PP8_plasma_abundance.csv")

# Read the csv file for the fetal sample of PP8 and calculate relative abundance
fetus8 <- read.csv("PP8_F_net2.csv")
sumareads <- sum(fetus8$Nreads)
fetus8$abundance <- (fetus8$`Nreads`/sumareads)* 100
rownames(fetus8) <- fetus8$X
fetus8$X <- NULL
fetus8$Rank <- NULL
fetus8 <- fetus8[order(fetus8$abundance, decreasing = TRUE), ]
write.csv2(fetus8,"PP8_fetus_abundance.csv")

# Read the csv file for the endometrial sample of PP8 and calculate relative abundance
endometri8 <-read.csv("PP8_C_net2.csv")
sumareads <- sum(endometri8$Nreads)
endometri8$abundance <- (endometri8$`Nreads`/sumareads)* 100
rownames(endometri8) <- endometri8$X
endometri8$X <- NULL
endometri8$Rank <- NULL
endometri8 <- endometri8[order(endometri8$abundance, decreasing = TRUE), ]
write.csv2(endometri8,"PP8_endometri_abundance.csv")

## --- FILTERING OF SAMPLES BY COMPARING CPMs BETWEEN SAMPLES AND CONTROLS

##Calculate CPM for extraction control and remove human reads
C_EXT <- C_EXT_S10_L001_fastp__kreport.txt %>%
 #divide reads by total reads and multiply by 1 million. Add the values to column called CPM
   mutate(CPM = (Nreads/1884401)*1000000)
#Remove Homo reads (human reads)
C_EXT <- C_EXT[C_EXT$GENUS != "Homo", ]
#Set genus names as rownames
rownames(C_EXT) <- C_EXT$GENUS
#Remove genus column
C_EXT$GENUS <- NULL
#Save filtered .csv
write.csv(C_EXT,"C_EXT_net2.csv")

##Calculate CPM for library control and remove human reads
C_LIB <- C_LIB_S5_L001_fastp__kreport.txt %>%
  mutate(CPM = (Nreads/22128)*1000000)
C_LIB <- C_LIB[C_LIB$GENUS != "Homo", ]
rownames(C_LIB) <- C_LIB$GENUS
C_LIB$GENUS <- NULL
#Save filtered .csv
write.csv(C_LIB,"C_LIB_net2.csv")


#Create a vector with all the genera present in any of the two negative controls
all_contam_genera <- union(rownames(C_LIB), rownames(C_EXT))

#For each genus present in the controls, take its CPM from C_EXT and C_LIB. 
#and calculate a reference CPM value
#Create a table with a column called genus from the vector all_contam_genera
contam_combined <- tibble(genus = all_contam_genera) %>%
  #Create new columns for cpm values of C_EXT and C_LIB
  mutate(
    cpm_contam1 = C_EXT[genus, 4],  # NA if the genus is not in contam1
    cpm_contam2 = C_LIB[genus, 4],  # NA if the genus is not in contam2
    # Reference value = max value between the two tables (ignoring NAs)
    cpm_ref     = pmax(cpm_contam1, cpm_contam2, na.rm = TRUE)
  ) %>%
  #Set genus as rownames
  column_to_rownames("genus")

#Save the reference contamination table
write.csv(contam_combined, "contam_combined.csv")

## --- Compare contam_combined table with sample tables (repeated as many times needed for each sample table)---

#Load the reference table (manually reviewed in excel) and the sample table
table1 <- read.csv("contam_combined_1_from_excel.csv")
table2 <- read.csv("PP7_V_net2.csv")

#Merge both tables according to genus.  
comparison <- merge(table1, table2, 
                    by.x = colnames(table1)[1], 
                    by.y = colnames(table2)[1], 
                    all.y = TRUE)  # keeps all the genera from table 2 even if not present in table 1


#Extract CPM values 
cpm_t1 <- comparison[, "cpm_ref"] # CPM reference column from table1
cpm_t2 <- comparison[, "CPM"]  # CPM column 5 from table2

#Mark as TRUE those genera where CPM table1 >= CPM table2
keep <- cpm_t1 <= cpm_t2

#Store the genera where sample CPM < reference CPM (to delete)
genus_delete <- comparison[!keep, 1]

#Eliminate those genera from table2
table2_filtered <- table2[!table2[, 1] %in% genus_delete, ]

#Save the filtered table
write.csv(table2_filtered, "PP7_V_filtered.csv", row.names = FALSE)

#Print a summary of the filtering results
cat("Originall rows in tabla2:", nrow(table2), "\n")
cat("Deleted rows:          ", nrow(table2) - nrow(table2_filtered), "\n")
cat("Conserved rows:         ", nrow(table2_filtered), "\n")

#Repeat this part for all the samples needed, changing the input file


##--- ALPHA DIVERSITY CALCULATION----

#Extract read counts from the vaginal sample
reads <- VAGINA$Nreads
#Calculate relative proportions (must add 1)
props <- reads / sum(reads)  # relative abundance

#Observed: total number of genera with at least 1 read
observed <- sum(reads > 0)
cat("Observed:", observed, "\n")

#Shannon index = -sum(pi * ln(pi))
#To measure richness of the sample
shannon <- -sum(props * log(props))
cat("Shannon:", shannon, "\n")

#Simpson index= 1-sum(pi^2)
#To measure dominance and diversity
simpson <- 1 - sum(props^2)
cat("Simpson:", simpson, "\n")

#Create a table with the three metrics
alpha_div <- data.frame(
  Metrics = c("Observed", "Shannon", "Simpson"),
  Value = c(observed, shannon, simpson)
)

