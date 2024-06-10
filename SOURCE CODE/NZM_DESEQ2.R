library(tximport)
library(DESeq2)
library(apeglm)
library(pryr)
library(tidyverse)
library(dplyr)

#DIFFERENTIAL EXPRESSION USING TEtranscripts data
data <- read.table("/Users/idakups/Desktop/Thesis/TEcountData/NZM/NZM.cntTable",header=T,row.names=1)
#groups <- factor(c(rep("CGroup",5),rep("TGroup",4))) #DOES THAT ASUME THE RIGHT ORDER
groups <- factor(c('CGroup','CGroup','CGroup','TGroup','TGroup','CGroup','CGroup','TGroup','TGroup'))
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


library(tximport)
#ZMIENILAM KOLEJNOSC W METADATA!!!
deseq_dir_path <- '/Users/idakups/Desktop/Thesis/RsemData/NZM'
deseq_samples <- list.files(deseq_dir_path)
deseq_metadata <- read.table(file.path(deseq_dir_path, 'metadata.txt'), header = TRUE)

deseq_sample_paths <- file.path(deseq_dir_path, deseq_samples[grep(".results$", deseq_samples)])
names(deseq_sample_paths) <- paste0("human", 1:9)
txi.rsem <- tximport(deseq_sample_paths, type = "rsem", txIn = FALSE, txOut = FALSE)

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

deseq_dds <- DESeqDataSetFromTximport(txi.rsem, colData = deseq_metadata, design = ~ condition  )
deseq_dds <- DESeq(deseq_dds)
deseq_res <- results(deseq_dds)
deseq_res

deseq_res

deseq_res <- na.omit(deseq_res)

logGreaterThan20 <- deseq_res[abs(deseq_res$log2FoldChange)>20,]


logGreaterThan20


#extract the differentially expressed ERV's
ERV_df <- deseq_res[str_starts(rownames(deseq_res), 'Hsap38'), ]
gene_df <- deseq_res[str_starts(rownames(deseq_res), 'ENS'), ]
ERV_df

#write the data to file
file_path = '/Users/idakups/Desktop/Thesis/RsemData/NZM/NZM_deseq_erv_full_results.tsv'
write.table(ERV_df, file = file_path, sep = "\t", row.names = TRUE)


ERV_diff_expr <- nrow(ERV_df[ERV_df$padj < 0.1 & ERV_df$log2FoldChange> 0.6,]) #29
ERV_non_diff_expr <- nrow(ERV_df[ERV_df$padj >= 0.1 ,])

gene_diff_expr <- nrow(gene_df[gene_df$padj < 0.1 & gene_df$log2FoldChange  > 0.6,]) #941
gene_non_diff_expr <- nrow(gene_df[gene_df$padj >= 0.1,])


ERV_diff_expr
ERV_non_diff_expr
gene_diff_expr
gene_non_diff_expr

gene_diff_expr <- gene_df[gene_df$padj < 0.1 & abs(gene_df$log2FoldChange)> 0.6,]
ERV_diff_expr <- ERV_df[ERV_df$padj < 0.1 & abs(ERV_df$log2FoldChange)> 0.6,]

gene_diff_expr <- gene_diff_expr[, "log2FoldChange", drop = FALSE]

ERV_diff_expr <- ERV_diff_expr[, "log2FoldChange", drop = FALSE]

write.table(gene_diff_expr, file = '/Users/idakups/Desktop/Thesis/RsemData/NZM/NZM_genes_diff_expr.tsv', sep='\t', col.names = NA)
write.table(ERV_diff_expr, file = '/Users/idakups/Desktop/Thesis/RsemData/NZM/NZM_ERV_diff_expr.tsv', sep='\t', col.names = NA)


data= matrix(c(ERV_diff_expr, gene_diff_expr, ERV_non_diff_expr, gene_non_diff_expr), ncol=2, nrow= 2, byrow=TRUE)

# specify the column names and row names of matrix
colnames(data) = c('ERV','non.ERV')
rownames(data) <- c('diff.expr','non.diff.expr')

# assign to table
final=as.table(data)


fisher.test(final)
chisq.test(final)

#COUNT HOW MANY UP AND HOW MANY DOWN REGULATED
ERV <- ERV_df[ERV_df$padj < 0.1 & abs(ERV_df$log2FoldChange)> 0.6,]
ERV_up <- nrow(ERV[ERV$log2FoldChange > 0,])
ERV_up 
ERV_down <- nrow(ERV[ERV$log2FoldChange < 0,])
ERV_down


gene <- gene_df[gene_df$padj < 0.1 & abs(gene_df$log2FoldChange)> 0.6,]
#VOLCANO PLOT
library(ggrepel)
deseq_res$diffexpressed <- 'NO'
deseq_res$diffexpressed[deseq_res$log2FoldChange > 0.6 & deseq_res$padj < 0.1] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
deseq_res$diffexpressed[deseq_res$log2FoldChange < -0.6 & deseq_res$padj < 0.05] <- "DOWN"

deseq_res$gene_id <- rownames(deseq_res)
deseq_res

deseq_res$delabel <- NA
deseq_res$delabel[deseq_res$diffexpressed != "NO"] <- deseq_res$gene_id[deseq_res$diffexpressed != "NO"]


ggplot(data=deseq_res, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  ggtitle('NZM') +
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
file_path <- "/Users/idakups/Desktop/Thesis/RsemData/NZM/ERV_names.txt"


writeLines(ERV_names, file_path)


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

