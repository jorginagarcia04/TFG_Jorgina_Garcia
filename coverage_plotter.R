## --- COVERAGE PLOT FOR POSITIVE SAMPLES - PAPER 1 --- ##
## IDENTIFYING CONTAMINANTS IN METAGENOMICS AT FAMILY/GENUS LEVEL ##
## R version 4.2.1 (2022-06-23)
## Copyright: Marta Ibañez Lligoña (marta.ibanez@vhir.org)
## Adapted by: Jorgina Garcia Larroy

## --- COVERAGE PLOT --- ##

##Install packages if needed
BiocManager::install(c("Rsamtools", "GenomicAlignments", "Gviz")) #(if not installed)

library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)
#library(ggalt)
#library("ggbreak")  
#library(ggpubr)
#library(ggrepel)
#library("hexbin")
#library(patchwork)

## Set working directory
setwd("C:/Users/jorgi/OneDrive - Universitat Pompeu Fabra/BIO HUMANA/4t BIO HUMANA/TFG/Analysis/Kraken")

## Define input BAM file
bam_file <- ("C:/Users/jorgi/OneDrive - Universitat Pompeu Fabra/BIO HUMANA/4t BIO HUMANA/TFG/Analysis/Kraken/PP8_pbivia_sorted.bam")


##Define reference genome
ref_PB <- "NZ_AP038912.1" #Prevotella bivia reference genome


# Function to generate coverage plots
# Calculates coverage per base and generates a plot for the selected chromosome
plot_genome_coverage <- function(bam_file, ref_genome, sample_label) {
 
  #Read mapped sequencing reads from the BAM file
  alignments <- readGAlignments(bam_file)
  #Calculate coverage across the reference sequence
  cov_data <- coverage(alignments)
  
  #Verify that the selected chromosome is present in the alignment file
  if (!ref_genome %in% names(cov_data)) {
    stop(paste("Region", ref_genome, "not found in coverage object."))
  }
  
  # Extract coverage for selected chromosome
  genome_cov <- as.numeric(cov_data[[ref_genome]])
  
  # Create a data frame containing genomic coordinates, corresponding to coverage values
  cov_df <- data.frame(
    position = seq_along(genome_cov),
    coverage = genome_cov
  )

  
  # Generate coverage plot
  p <- ggplot(cov_df, aes(x = position, y = coverage)) +
    geom_line(color = "steelblue3") + 
    #Define plot totle and axis labels
    labs(
      title = paste(sample_label),
      x = "Genomic Position",
      y = "Read depth"
    ) +
    #Display x-axis labels every 300 kbases
    scale_x_continuous(breaks = seq(0, max(cov_df$position), by = 300000)) + 
    #display coverage to 250x for better visualisation
    theme_bw() + ylim(0,250) +
    theme(
      #adjust text size
      axis.text.x = element_text(size = 10),
      axis.title = element_text(size = 12),
      title = element_text(size = 12)
    )
  
  # Return both plot and coverage data
  return(list(
    plot = p,
    coverage_data = cov_df,
    mean_coverage = mean(genome_cov)
  ))
}
#Generate coverage plot for Prevotella bivia
prevotella <- plot_genome_coverage(bam_file, ref_PB, "PP8_PT")

#Display coverage plot
prevotella$plot
#Calculate summary coverage statistics
prevotella.mean_coverage <- data.frame(
  sample = "PP8_PT_PB",
  #mean sequencing depth across the chromosome
  mean_coverage = prevotella$mean_coverage,
  #Breadth of ceverage calculated: percentage of genomic positions
  #covered by at least one read
  percent_covered = (sum(prevotella$coverage_data$coverage > 0) / nrow(prevotella$coverage_data)) * 100
)

