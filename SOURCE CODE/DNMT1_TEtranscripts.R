library(tidyverse)

data <- read.table("/Users/idakups/Desktop/Thesis/TEcountData/DNMT1/DNMT1.cntTable",header=T,row.names=1)

groups <- factor(c(rep("CGroup",2),rep("TGroup",2)))
min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds$groups = relevel(dds$groups,ref="CGroup")
dds <- DESeq(dds)
res <- results(dds)

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="groups_TGroup_vs_CGroup", type="apeglm")


plotMA(res)
plotMA(resLFC)
sum(res$padj < 0.1, na.rm=TRUE) #334
sum(res$padj < 0.05, na.rm=TRUE) #235


#DESEQ 
library(tximport)
library(DESeq2)
library(apeglm)
library(pryr)

dir_path <- '/Users/idakups/Desktop/Thesis/RsemData/DNMT1'
samples <- list.files(dir_path)
metadata <- read.table(file.path(dir_path, 'DNMT1_metadata.txt'), header = TRUE)

sample_paths <- file.path(dir_path, samples[grep(".results$", samples)])
names(sample_paths) <- paste0("mouse", 1:4)
txi.rsem <- tximport(sample_paths, type = "rsem", txIn = FALSE, txOut = FALSE)

#filter out zero length transcripts
txi.rsem$abundance <-
  txi.rsem$abundance[apply(txi.rsem$length,
                           1,
                           function(row) all(row !=0 )),]

txi.rsem$counts <-
  txi.rsem$counts[apply(txi.rsem$length,
                        1,
                        function(row) all(row !=0 )),]

txi.rsem$length <-
  txi.rsem$length[apply(txi.rsem$length,
                        1,
                        function(row) all(row !=0 )),]

# create an object storing the read counts and intermediate estimated quantities
dds2 <- DESeqDataSetFromTximport(txi.rsem, colData = metadata, design = ~ condition  )

dds2 <- DESeq(dds)
res2 <- results(dds)
res2

res2 <- na.omit(res)
res2

#PLOT LOG2FOLDCHANGES RSEM VS TETRANSCRIPT
# Get the intersection of row names
common_row_names <- intersect(rownames(res), rownames(res2))
# Subset both data frames to keep only the rows with common row names
res_common <- res[common_row_names, ]
deseq_res_common <- res2[common_row_names, ]

names(deseq_res_common)[1] <- 'baseMean2'
names(deseq_res_common)[2] <- 'log2FoldChange2'
names(deseq_res_common)[3] <- 'lfcSE2'
names(deseq_res_common)[4] <- 'stat2'
names(deseq_res_common)[5] <- 'pvalue2'
names(deseq_res_common)[6] <- 'padj2'

res_common <- as.data.frame(res_common)
deseq_res_common <- as.data.frame(deseq_res_common)

merged_df <- merge(res_common, deseq_res_common, 
                   by = 'row.names', all = TRUE)
#res_common x, dnmt_res_common y
g <- ggplot(merged_df, aes(x = log2FoldChange, y= log2FoldChange2) ) + geom_point() + xlim(-12, 10) + geom_abline(intercept = 0, slope = 1, color = "red")
g + labs(title="Log2FoldChange RSEM vs TEtranscript", y="RSEM Log2FoldChange", x="TEtranscript Log2FoldChange")

library(ggplot2)
library(ggExtra)

# Create the ggplot object without marginal distributions
g <- ggplot(merged_df, aes(x = log2FoldChange, y = log2FoldChange2)) +
  geom_point() +
  xlim(-12, 10) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = "Log2FoldChange RSEM vs TEtranscript", y = "RSEM Log2FoldChange", x = "TEtranscript Log2FoldChange") 
p <- ggMarginal(g, type = 'histogram')
p

p2  <- ggMarginal(g, type = 'density')
p2



#GET OVEREXPRESSED ERVS THRESHOLD 
threshold <- 0.05
res3 <- res[grep("^Hsap38", rownames(res)), ]
res3

overexpressed_ERVs <- res3[res3$pvalue <= threshold, ]
overexpressed_ERVs

ERV_names <- as.character(rownames(overexpressed_ERVs))
ERV_names

threshold <- 0.1

res3 <- na.omit(res3)
overexpressed_ERVs2 <- res3[res3$padj <= threshold, ]
overexpressed_ERVs2  <- overexpressed_ERVs2[abs(overexpressed_ERVs2$log2FoldChange) > 0,]
overexpressed_ERVs2

ERV_names2 <- as.character(rownames(overexpressed_ERVs2))
ERV_names2

# Specify the file path where you want to save the text file
file_path <- "/Users/idakups/Desktop/Thesis/RsemData/DNMT1/ERV_namesTEtranscript.txt"

# Write the ERV_names vector to the text file
writeLines(ERV_names2, file_path)
