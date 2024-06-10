library(tximport)
library(DESeq2)
library(apeglm)
library(pryr)
library(tidyverse)
library(dplyr)

#TETRANSCRIPT PART 
data <- read.table("DNMT1.cntTable",header=T,row.names=1)
groups <- factor(c(rep("TGroup",2),rep("CGroup",2)))
min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds$groups = relevel(dds$groups,ref="CGroup")
dds <- DESeq(dds)
res <- results(dds)
write.table(res, file="DNMT1_gene_TE_analysis.txt", sep="\t",quote=F)
resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) &         (abs(res$log2FoldChange)> 0.000000)), ]
write.table(resSig, file="DNMT1_sigdiff_gene_TE.txt",sep="\t", quote=F)







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
dds <- DESeqDataSetFromTximport(txi.rsem, colData = metadata, design = ~ condition  )

dds <- DESeq(dds)
res <- results(dds)
res 

res <- na.omit(res)
res #gene expression with pvalue 

pval_threshold <- 0.1
# Assuming 'res' is a DESeqResults object
# Extract the results and convert it to a dataframe


# Filter rows based on row names starting with "ENS"
gene_df <- res[str_starts(rownames(res), 'ENS'), ]
ERV_df <- res[str_starts(rownames(res), 'Hsap38'), ]

# write ERV deseq results
file_path = '/Users/idakups/Desktop/Thesis/RsemData/DNMT1_deseq_erv_full_results.tsv'
write.table(ERV_df, file = file_path, sep = "\t", row.names = TRUE)

gene_diff_expr <- nrow(gene_df[gene_df$padj < 0.1 & gene_df$log2FoldChange < -0.6,]) #3405 up 2008 down
gene_non_diff_expr <- nrow(gene_df[gene_df$padj >= 0.1,])

ERV_diff_expr <- nrow(ERV_df[ERV_df$padj < 0.1 & ERV_df$log2FoldChange < -0.6,]) #102 up 34 down 
ERV_non_diff_expr <- nrow(ERV_df[ERV_df$padj >= 0.1 ,])

gene_diff_expr
gene_non_diff_expr
ERV_diff_expr 
ERV_non_diff_expr

gene_diff_expr <- gene_df[gene_df$padj < 0.1 & abs(gene_df$log2FoldChange)> 0.6,]
ERV_diff_expr <- ERV_df[ERV_df$padj < 0.1 & abs(ERV_df$log2FoldChange)> 0.6,]

gene_diff_expr <- gene_diff_expr[, "log2FoldChange", drop = FALSE]

ERV_diff_expr <- ERV_diff_expr[, "log2FoldChange", drop = FALSE]

write.table(gene_diff_expr, file = '/Users/idakups/Desktop/Thesis/RsemData/DNMT1/DNMT1_genes_diff_expr.tsv', sep='\t', col.names = NA)
write.table(ERV_diff_expr, file = '/Users/idakups/Desktop/Thesis/RsemData/DNMT1/DNMT1_ERV_diff_expr.tsv', sep='\t', col.names = NA)

# create matrix with 4 columns and 4 rows
data= matrix(c(ERV_diff_expr, gene_diff_expr, ERV_non_diff_expr, gene_non_diff_expr), ncol=2, nrow= 2, byrow=TRUE)

# specify the column names and row names of matrix
colnames(data) = c('ERV','non.ERV')
rownames(data) <- c('diff.expr','non.diff.expr')

# assign to table
final=as.table(data)

mosaicplot(final,
           main = "Mosaic plot",
           color = TRUE
)

# display
final


fisher.test(final)
chisq.test(final)

#MAKE VOLCANO PLOT
library(ggrepel)
res$diffexpressed <- 'NO'
res$diffexpressed[res$log2FoldChange > 0.6 & res$padj < 0.1] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res$diffexpressed[res$log2FoldChange < -0.6 & res$padj < 0.05] <- "DOWN"

res$gene_id <- rownames(res)
res

res$delabel <- NA
res$delabel[res$diffexpressed != "NO"] <- res$gene_id[res$diffexpressed != "NO"]


ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  ggtitle('DNMT1+3B') +
  theme(plot.title = element_text(hjust=0.5, size = 20, face = 'bold')) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

plotMA(res)

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_B_vs_A", type="apeglm")
resLFC

plotMA(resLFC)


txi.rsem2 <-txi.rsem

txi.rsem2$abundance <- txi.rsem2$abundance[grep("^Hsap38", rownames(txi.rsem2$abundance)), ]
txi.rsem2$counts <- txi.rsem2$counts[grep("^Hsap38", rownames(txi.rsem2$counts)), ]
txi.rsem2$length <- txi.rsem2$length[grep("^Hsap38", rownames(txi.rsem2$length)), ]

dds2 <- DESeqDataSetFromTximport(txi.rsem2, colData = metadata, design = ~ condition  )

dds2 <- DESeq(dds2)
res2 <- results(dds2)
res2

plotMA(res2)

res2 <- na.omit(res2)

#MAKE VOLCANOPLOT SHOWING OVEREXPRESSED ERVS
res2$differentially.expressed <- 'NO'
res2$differentially.expressed[res2$log2FoldChange > 0.6 & res2$padj < 0.1] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res2$differentially.expressed[res2$log2FoldChange < -0.6 & res2$padj < 0.05] <- "DOWN"

res2$gene_id <- rownames(res2)
res2

res2$delabel <- NA
res2$delabel[res2$differentially.expressed != "NO"] <- res2$gene_id[res2$differentially.expressed != "NO"]


ggplot(data=res2, aes(x=log2FoldChange, y=-log10(padj), col=differentially.expressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size = 3) +
  ggtitle('DNMT1+3B - ERV expression') +
  theme(plot.title = element_text(hjust=0.5, size = 20, face = 'bold')) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red")


ggplot(data=res2, aes(x=log2FoldChange, y=-log10(padj), col=differentially.expressed)) +
  geom_point() + 
  theme_minimal() +
  ggtitle('DNMT1+3B - ERV expression') +
  theme(plot.title = element_text(hjust=0.5, size = 20, face = 'bold')) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red")


#PRINT OVEREXPRESSED ERV's 
ERV_diff_expr <- ERV_df[ERV_df$padj < 0.1 & abs(ERV_df$log2FoldChange)> 0.6,]
gene_diff_expr <- gene_df[gene_df$padj < 0.1 & abs(gene_df$log2FoldChange)> 0.6,]


#overexpressed_ERVs <- res2[res2$padj < 0.1, ]
#overexpressed_ERVs <- overexpressed_ERVs[abs(overexpressed_ERVs$log2FoldChange) > 0.6,]
#overexpressed_ERVs

#r1 <- rownames(ERV_diff_expr)
#r2 <- rownames(overexpressed_ERVs)

#intersect(r1, r2)
gene_names <- as.character(rownames(gene_diff_expr))

ERV_names <- as.character(rownames(ERV_diff_expr))
ERV_names

# Specify the file path where you want to save the text file
file_path <- "/Users/idakups/Desktop/Thesis/RsemData/DNMT1/ERV_names.txt"
file_path_gene <- "/Users/idakups/Desktop/Thesis/RsemData/DNMT1/gene_names.txt"

# Write the ERV_names vector to the text file
writeLines(ERV_names, file_path)
writeLines(gene_names, file_path_gene)


#PLOT LOG2FOLDCHANGES RSEM VS TETRANSCRIPT
# Get the intersection of row names
common_row_names <- intersect(rownames(res), rownames(deseq_res))

# Subset both data frames to keep only the rows with common row names
res_common <- res[common_row_names, ]
deseq_res_common <- deseq_res[common_row_names, ]

logGreaterThan20 <- deseq_res_common[abs(deseq_res_common$log2FoldChange)>20,]
logGreaterThan20



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

