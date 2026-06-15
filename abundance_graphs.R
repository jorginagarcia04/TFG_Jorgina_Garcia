## --- ABUNDANCE GRAPH CREATION FROM SAMPLES --- ##
## Jorgina Garcia Larroy

## --- Set working directory --- #
setwd("C:/Users/jorgi/OneDrive - Universitat Pompeu Fabra/BIO HUMANA/4t BIO HUMANA/TFG/Analysis/Kraken/kraken_reads")

## --- Load libraries, if necessary --- ##
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(RColorBrewer)

## ---------------------------------- PP7 SAMPLES ---------------------------------- ##
## --- ABUNDANCE GRAPH VAGINA --- ##
## --- Load data --- ##
vagina <- read.csv("PP7_V_net2.csv")

## Compute relative abundance (not taking into account unclassified or human reads)
sumareads <- sum(vagina$Nreads)
vagina$abundance <- (vagina$`Nreads`/sumareads)* 100
#Filter genera by abundance. Keeping abundance higher than 0,3 to reduce noise
vagina2 <- vagina[vagina$abundance > 0.3,]

## The genus with an abundance lower than 0,3 are grouped and the total abundance is added
other_abundance <- sum(vagina$abundance[vagina$abundance <= 0.3])

#Keep only genus name and abundaance columns
vagina3 <- vagina2[, c("X", "abundance")]
#Create new row for the "Other category" with the addition of abundances
other <- data.frame(
  X = "Other",
  abundance = other_abundance
)

#Combine the filtered table with the "Other" row
vagina_plot <- rbind(vagina3, other)

#Create a barplot of relative abuundances ordered from highest to lowest
ggplot(vagina_plot, aes(x = reorder(X, -abundance), y = abundance, fill = X)) + geom_col() + theme_bw() + scale_x_discrete() + xlab("Genus") + ylab("Relative abundance (%)") +
  #Set y axis from 0 to 100 with breaks every 5
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=10), axis.title = element_text(size=13), legend.position = "none") + scale_y_continuous( limits = c(0,100),
    breaks = seq(0, 100, by = 5)
  ) +
  #assign colours to each genus
  scale_fill_manual(values = c("#79CDCD", "#EEE9BF", "#98AFC7", "#B4EEB4","#FFC1C1", "#4F94CD", "#C6E2FF","#CD6889","#6698FF", "#FF6A6A", "pink", "#698B69", "#FFA07A")) + 
  scale_color_manual(values = c("#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF", "#FFFFFF"))

#Save the figure as .png
ggsave("Figure1_vaginal_microbiome.png",
       width = 8,
       height = 6,
       dpi = 600)


## --- HEATMAP PP7----

#Load filtered tables of the remaining PP7 samples
placenta <- read.csv("PP7_PT_net2.csv")
plasma <- read.csv("PP7_PM_net2.csv")
endometri <- read.csv("PP7_E_net2.csv")

#Calculate relative abundances and filter genera >0,15
sumareads <- sum(placenta$Nreads)
placenta$abundance <- (placenta$`Nreads`/sumareads)* 100
placenta7 <- placenta[placenta$abundance > 0.15,]

#Filter vagina genera >0,15
vagina7 <- vagina[vagina$abundance > 0.15,]

#Calculate relative abundance and filter the genera >0,15%
sumareads <- sum(plasma$Nreads)
plasma$abundance <- (plasma$`Nreads`/sumareads)* 100
plasma7 <- plasma[plasma$abundance > 0.15,]

#Calculate relative abundance and filter the genera >0,15%
sumareads <- sum(endometri$Nreads)
endometri$abundance <- (endometri$`Nreads`/sumareads)* 100
endometri7 <- endometri[endometri$abundance > 0.15,]

#Keep only the genus name (column X) and abundance for each sample
placenta7 <- placenta7[, c("X", "abundance")]
vagina7 <- vagina7[, c("X", "abundance")]
plasma7 <- plasma7[, c("X", "abundance")]
endometri7 <- endometri7[, c("X", "abundance")]

#Rename abundance columns with sample names
colnames(vagina7)[colnames(vagina7) == "abundance"] <- "PP7_V"
colnames(placenta7)[colnames(placenta7) == "abundance"] <- "PP7_PT"
colnames(endometri7)[colnames(endometri7) == "abundance"] <- "PP7_E"
colnames(plasma7)[colnames(plasma7) == "abundance"] <- "PP7_PM"

#Merge all PP7 tables by genus name, keeping all the genera
merge1 <- merge(vagina7, placenta7, by = "X", all = TRUE)
merge2 <- merge(merge1, endometri7, by = "X", all = TRUE)
merge3 <- merge(merge2, plasma7, by = "X", all = TRUE)

#Replace NAs with 0
merge3[is.na(merge3)] <- 0
#Calculate cumulative abundance across all samples for each genus
merge3$Total <- rowSums(merge3[, -1])
#Order genera by decreasing abundance
merge3 <- merge3[order(merge3$Total, decreasing = TRUE), ]
#Keep the genera with cumulative abundance >1,5%
merge4 <- merge3[merge3$Total > 1.5,]
#Set genus names as row names
rownames(merge4) <- merge4$X
merge4$X <- NULL

#Remove the "Total" column before plotting
heatmap_abundance <-merge4[,-5]

library(pheatmap)
#Log10 transformation of the abundance values + 0.01 to avoid log(0)
heatmap_log <- log10(heatmap_abundance + 0.01)

# Define genera classified as low confidence 
low_confidence <- c("Burkholderia", "Escherichia", "Paraburkholderia", "Acinetobacter", "Cutibacterium",
                    "Klebsiella", "Rhodococcus", "Streptomyces", "Pseudomonas", "Cupriavidus", "Bdellovibrio",
                    "Corynebacterium", "Pendulispora", "Salmonella", "Malassezia",
                    "Polynucleobacter", "Yersinia", "Nocardioides", "Mycoavidus", "Mycobacterium", 
                    "Streptococcus", "Enterococcus", "Flavobacterium")

# Create dataframe with row annotation containing confidence level
row_annot <- data.frame(
  Confidence = ifelse(rownames(heatmap_log) %in% low_confidence,
                      "Low confidence",
                      "High confidence")
)
rownames(row_annot) <- rownames(heatmap_log)

# Define colors for the annoted confidence level
annot_colors <- list(
  Confidence = c("Low confidence" = "grey80",
                 "High confidence" = "#FFA07A")
)

# Plot the heatmap with confidence annotation
pheatmap(
  heatmap_log,
  angle_col = 45, #rotate column labels 45º
  border_color = NA, #eliminate border color
  cellwidth = 40,
  legend_breaks = c(-2, -1, 0, 1, 1.5, 1.85),
  #legend labels show corresponding percentage values for relative abundance
  legend_labels = c("0%", "0.1%", "1%", "10%", "30%", "70%"),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = colorRampPalette(c("white", "lightblue", "darkblue"))(100),
  fontsize_row = 9,
  fontsize_col = 11,
  annotation_row = row_annot,          # add lateral confidence bar
  annotation_colors = annot_colors      # add colors for the annotation bar
)

## ---------------------------------- PP8 SAMPLES ---------------------------------- ##
## --- HEATMAP PP8--- ##

#Load data for PP8 samples 
placenta8 <-read.csv("PP8_PT_net2.csv")
plasma8 <- read.csv("PP8_PM_net2.csv")
fetus8 <- read.csv("PP8_F_net2.csv")
endometri8 <-read.csv("PP8_C_net2.csv")

#Calculate relative abundances and filter genera >0,2% relative abundance
sumareads <- sum(placenta8$Nreads)
placenta8$abundance <- (placenta8$`Nreads`/sumareads)* 100
placenta8 <- placenta8[placenta8$abundance > 0.2,]

sumareads <- sum(plasma8$Nreads)
plasma8$abundance <- (plasma8$`Nreads`/sumareads)* 100
plasma8 <- plasma8[plasma8$abundance > 0.2,]

sumareads <- sum(fetus8$Nreads)
fetus8$abundance <- (fetus8$`Nreads`/sumareads)* 100
fetus8 <- fetus8[fetus8$abundance > 0.2,]

sumareads <- sum(endometri8$Nreads)
endometri8$abundance <- (endometri8$`Nreads`/sumareads)* 100
endometri8 <- endometri8[endometri8$abundance > 0.2,]

#Keep only genus name and abundance columns for each sample
placenta8 <- placenta8[, c("X", "abundance")]
fetus8 <- fetus8[, c("X", "abundance")]
plasma8 <- plasma8[, c("X", "abundance")]
endometri8 <- endometri8[, c("X", "abundance")]

#Rename columns with sample names
colnames(fetus8)[colnames(fetus8) == "abundance"] <- "PP8_F"
colnames(placenta8)[colnames(placenta8) == "abundance"] <- "PP8_PT"
colnames(endometri8)[colnames(endometri8) == "abundance"] <- "PP8_C"
colnames(plasma8)[colnames(plasma8) == "abundance"] <- "PP8_PM"

#Merge all PP8 tables by genus name, keeping all genera
merge81 <- merge(fetus8, placenta8, by = "X", all = TRUE)
merge82 <- merge(merge81, endometri8, by = "X", all = TRUE)
merge83 <- merge(merge82, plasma8, by = "X", all = TRUE)

#Replace NAs with 0s
merge83[is.na(merge83)] <- 0
#Calculate cumulative abundance and order by decreasing total abundance
merge83$Total <- rowSums(merge83[, -1])
merge83 <- merge83[order(merge83$Total, decreasing = TRUE), ]

#Keep only genera with cumulative abundance >1,5%
merge84 <- merge83[merge83$Total > 1,5,]
#Set genus as rownames and eliminate the genus column (called X)
rownames(merge84) <- merge84$X
merge84$X <- NULL

#remove the "Total" column before plotting
heatmap_abundance2 <-merge84[,-5]
#Log10 transformation of the abundance values
heatmap_log2 <- log10(heatmap_abundance2 + 0.01)


library(pheatmap)

#Define confidence levels for the genera of the plot
high_confidence <- c("Prevotella", "Klebsiella", "Lactobacillus", "Veillonella", "Gardnerella")   # aparecen en 1 muestra
med_confidence <- c("Saalmonella", "Enterobacter", "Lawsonella")    # aparecen en 2-3 muestras
low_confidence <- c("Burkholderia", "Escherichia", "Paraburkholderia", "Acinetobacter", "Cutibacterium",
                    "Rhodococcus", "Streptomyces", "Pseudomonas", "Cupriavidus", "Bdellovibrio",
                    "Corynebacterium", "Pendulispora", "Vibrio", "Malassezia",
                    "Polynucleobacter", "Yersinia", "Nocardioides", "Mycoavidus", "Mycobacterium", 
                    "Streptococcus", "Enterococcus", "Flavobacterium", "Staphylococcus",
                    "Enterococcus", "Rhodococcus", "Micrococcus","Kocuria", "Moraxella",
                    "Aquirufa")    

# Crea row annotation dataframe with the confidence level annotation
row_annot2 <- data.frame(
  Confidence = case_when(
    rownames(heatmap_log2) %in% low_confidence ~ "Low confidence",
    rownames(heatmap_log2) %in% med_confidence ~ "Medium confidence",
    TRUE ~ "High confidence" #all remaining genera are classified as High confidence
  )
)
rownames(row_annot2) <- rownames(heatmap_log2)

# Define colors for the three confidence levels
annot_colors <- list(
  Confidence = c("Low confidence"    = "grey80",
                 "Medium confidence" = "#FFE7BA",
                 "High confidence"   = "#FFA07A")
)
# Pheatmap plot with annotation
pheatmap(
  heatmap_log2,
  angle_col = 45,
  border_color = NA,
  cellwidth = 40,
  #custom breaks for PP8
  legend_breaks = c(-2.0 , -0.68 , 0.00,  0.70 , 1.00 ,  1.54),
  legend_labels = c("0%", "0.2%", "1%", "5%", "10%", "35%"),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = colorRampPalette(c("white", "lightblue", "darkblue"))(100),
  fontsize_row = 9,
  fontsize_col = 11,
  annotation_row = row_annot2,          # adds lateral bar of annotations
  annotation_colors = annot_colors      # adds colours
)


## --- JOINT HEATMAP PP7 AND PP8 ---

#  Merge all comparable samples from PP7 and PP8 by genus name
heatmap_conjunt <- endometri7 %>%
  full_join(plasma7, by = "X") %>%
  full_join(placenta7, by = "X") %>%
  full_join(endometri8,  by = "X") %>%
  full_join(plasma8, by = "X") %>%
  full_join(placenta8, by = "X")

#Replace NAs with 0
heatmap_conjunt[is.na(heatmap_conjunt)] <- 0
#Set genus as rownames
rownames(heatmap_conjunt) <- heatmap_conjunt$X
heatmap_conjunt$X <- NULL
#Calculate cumulative abundance across all 6 samples and order by decreasing total abundance
heatmap_conjunt$Total <- rowSums(heatmap_conjunt[, -1])
heatmap_conjunt <- heatmap_conjunt[order(heatmap_conjunt$Total, decreasing = TRUE), ]
#Keep only genera with cumulative abundance > 2,5%
heatmap_conjunt <- heatmap_conjunt[heatmap_conjunt$Total > 2.5,]

#Remove the "Total" column
heatmap_conjunt2 <-heatmap_conjunt[,-7]

# Define low confidence genera
low_confidence <- c("Burkholderia", "Escherichia", "Paraburkholderia", "Acinetobacter", "Cutibacterium",
                   "Rhodococcus", "Streptomyces", "Pseudomonas", "Cupriavidus", "Bdellovibrio",
                   "Corynebacterium", "Pendulispora", "Vibrio", "Malassezia",
                   "Polynucleobacter", "Yersinia", "Nocardioides", "Mycoavidus", "Mycobacterium", 
                   "Streptococcus", "Enterococcus", "Flavobacterium", "Staphylococcus",
                   "Enterococcus", "Rhodococcus", "Micrococcus","Kocuria", "Moraxella",
                   "Aquirufa", "Enterobacter", "Salmonella")

# Create row annotation dataframe with two confidence levels
row_annot <- data.frame(
  Confidence = ifelse(rownames(heatmap_conjunt2) %in% low_confidence, 
                      "Low confidence",
                      "High confidence")
)
rownames(row_annot) <- rownames(heatmap_conjunt2)

# Create column annotation dataframe indicating patient PP7 or PP8
col_annot <- data.frame(
  Patient = c("PP7", "PP7", "PP7", "PP8", "PP8", "PP8")
)
rownames(col_annot) <- colnames(heatmap_conjunt2)

# Define colors for the annotations
annot_colors <- list(
  Patient = c("PP7" = "#B4EEB4",
              "PP8" = "#66CDAA"),
  Confidence = c("Low confidence"  = "grey80",
                 "High confidence" = "#FFA07A")
)
#Log10 transformation of abundance values
heatmap_conjunt_log <- log10(heatmap_conjunt2 + 0.01)

# Heatmap with confidence and patient annotations
pheatmap(
  heatmap_conjunt_log,
  angle_col = 45,
  border_color = NA,
  cellwidth = 40,
  legend_breaks = c(-2, -1, 0, 1, 1.5, 1.85),
  legend_labels = c("0%", "0.1%", "1%", "10%", "30%", "70%"),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = colorRampPalette(c("white", "lightblue", "darkblue"))(100),
  fontsize_row = 9,
  fontsize_col = 11,
  annotation_row = row_annot,#lateral confidence bar
  annotation_col = col_annot, #top patient bar
  annotation_colors = annot_colors
)

