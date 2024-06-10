library(methylKit)
library(tidyverse)
library(genomation)

dir_path <- "/Users/idakups/Desktop/Thesis/methylation/TET2"
samples <- list.files(dir_path)
sample_cov_paths <- file.path(dir_path, samples[grep("_cov.txt$", samples)])
sample_cov_paths_list <- as.list(sample_cov_paths)

obj_cov <- methRead(location = sample_cov_paths_list, 
                    sample.id = list('cont1', 'cont2', 'mut1', 'mut2', 'mut3'),
                    assembly = 'mm10', 
                    treatment = c(0, 0, 1, 1, 1),
                    pipeline = "bismarkCoverage")


sample.id = list('cont1', 'cont2', 'mut1', 'mut2', 'mut3')

sample_data_list = list()
col_names <- c('chr', 'start', 'end', 'meth_percentage', 'numCs', 'numTs')

for (i in seq_along(sample.id)) {
  print(sample.id[i])
  print(sample_cov_paths_list)
  sample_table <- read.table(sample_cov_paths_list[[i]], col.names = col_names)  # Extract file path using double brackets
  sample_data_list[[sample.id[[i]]]] <- sample_table  # Use double brackets to assign values by name
}


sample_data_list[['cont1']] #if you wantt o access the first element of the list you have to use double brackets

sample_data_list[[1]]

fut8 <- c(12, 77237400, 77240801)
gata2_v1 <- c(6,88199000,88200600)
gata2_v2 <- c(6, 88193600, 88195001)
sirt5 <- c(6,13573200,13576001)



for (gene_name in names(filter_conditions)) {
  filter_cond <- filter_conditions[[gene_name]]
  sample_table <- subset(sample_table, chr == filter_cond[1] & start >= filter_cond[2] & end <= filter_cond[3])
  sample_data_list[[paste0(sample.id[[i]], "_", gene_name)]] <- sample_table
}
}
