library(GenomicFeatures)
library(tidyverse)

chrom_info <- getChromInfoFromBiomart(dataset = 'mmusculus_gene_ensembl',
                        host = "https://nov2020.archive.ensembl.org") %>%
  dplyr::arrange(chrom) 

chrom_info

# Specify the file path where you want to save the output
file_path <- "/Users/idakups/Desktop/Thesis/bedGraphtoBigWig/chromosome.txt"

# Write the data to a text file
write.table(chrom_info, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

# Print a message indicating successful file saving
cat("Chromosome information saved to:", file_path, "\n")
