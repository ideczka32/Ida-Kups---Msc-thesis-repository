library(tximport)
library(DESeq2)
library(apeglm)
library(pryr)
library(tidyverse)
library(dplyr)

dir_path <- '/Users/idakups/Desktop/Thesis/RsemData/TET2'
samples <- list.files(dir_path)
metadata <- read.table(file.path(dir_path, 'TET2_metadata.txt'), header = TRUE)


sample_paths <- file.path(dir_path, samples[grep(".results$", samples)])
names(sample_paths) <- paste0("mouse", 1:5)
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
res <- na.omit(res)
res #gene expression with pvalue 

padj_threshold <- 0.1

# Filter rows based on row names starting with "ENS"
gene_df <- res[str_starts(rownames(res), 'ENS'), ]
ERV_df <- res[str_starts(rownames(res), 'Mmus38'), ]

#write the data to file
file_path = '/Users/idakups/Desktop/Thesis/RsemData/TET2/TET2_deseq_erv_full_results.tsv'
write.table(ERV_df, file = file_path, sep = "\t", row.names = TRUE)


gene_diff_expr <- nrow(gene_df[gene_df$padj < 0.1 & abs(gene_df$log2FoldChange) > 0.6,]) 
gene_non_diff_expr <- nrow(gene_df[gene_df$padj >= 0.1,])

ERV_diff_expr <- nrow(ERV_df[ERV_df$padj < 0.1 & abs(ERV_df$log2FoldChange)> 0.6,])
ERV_non_diff_expr <- nrow(ERV_df[ERV_df$padj >= 0.1 ,])


gene_diff_expr
gene_non_diff_exprśś
ERV_diff_expr
ERV_non_diff_expr

# create matrix with 4 columns and 4 rows
data= matrix(c(ERV_diff_expr, gene_diff_expr, ERV_non_diff_expr, gene_non_diff_expr), ncol=2, nrow= 2, byrow=TRUE)

# specify the column names and row names of matrix
colnames(data) = c('ERV','non.ERV')
rownames(data) <- c('diff.expr','non.diff.expr')

# assign to table
final=as.table(data)

# display
final


fisher.test(final) #fisher's test
chisq.test(final) #large sample sizes


#VOLCANO PLOT
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
  ggtitle('TET2') +
  theme(plot.title = element_text(hjust=0.5, size = 20, face = 'bold')) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

#SAVE ERVS 
ERV_diff_expr <- ERV_df[ERV_df$padj < 0.1 & abs(ERV_df$log2FoldChange)> 0.6,]
gene_diff_expr <- gene_df[gene_df$padj < 0.1 & abs(gene_df$log2FoldChange)> 0.6,]

ERV_names <- as.character(rownames(ERV_diff_expr))
ERV_names

# Specify the file path where you want to save the text file
file_path <- "/Users/idakups/Desktop/Thesis/RsemData/TET2/ERV_names.txt"


writeLines(ERV_names, file_path)

ERV_diff_expr <- nrow(ERV_df[ERV_df$padj < 0.1 & ERV_df$log2FoldChange> 0.6,])
ERV_diff_expr

gene_diff_expr <- nrow(gene_df[gene_df$padj < 0.1 & gene_df$log2FoldChange> 0.6,])
gene_diff_expr
