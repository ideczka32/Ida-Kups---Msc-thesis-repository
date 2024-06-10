library(tximport)
library(DESeq2)
library(apeglm)
library(pryr)

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


#only use ERV rows
txi.rsem2 <-txi.rsem

txi.rsem2$abundance <- txi.rsem2$abundance[grep("^Mmus38", rownames(txi.rsem2$abundance)), ]
txi.rsem2$counts <- txi.rsem2$counts[grep("^Mmus38", rownames(txi.rsem2$counts)), ]
txi.rsem2$length <- txi.rsem2$length[grep("^Mmus38", rownames(txi.rsem2$length)), ]

print(txi.rsem$counts)
print(dim(txi.rsem2$counts))

# create an object storing the read counts and intermediate estimated quantities
dds <- DESeqDataSetFromTximport(txi.rsem, colData = metadata, design = ~ condition  )
#dds$condition <- factor(dds$condition, levels = c("W", "M"))

dds2 <- DESeqDataSetFromTximport(txi.rsem2, colData = metadata, design = ~ condition  )
#dds2$condition <- factor(dds$condition, levels = c("W", "M"))

#pre-filtering - should we apply?
#smallestGroupSize <- 3 #we want to keep only if at leas 3 samples have expression greater than 10
#keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
#dds <- dds[keep,]

#STANDARD DIFFERENTIAL ANALYSIS
dds <- DESeq(dds)
res <- results(dds)
res

dds2 <- DESeq(dds2)
res2 <- results(dds2)
res2

#filter out rows with no expression
res <- na.omit(res)
res

res2 <- na.omit(res2)
res2

#logfold change shrinkage
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_W_vs_M", type="apeglm")
resLFC

resultsNames(dds2)
resLFC2 <- lfcShrink(dds2, coef="condition_W_vs_M", type="apeglm")
resLFC2

#How many adjusted p-values were lower than 0.1
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.05, na.rm=TRUE)

sum(res2$padj < 0.1, na.rm=TRUE) #5081 total
sum(res2$padj < 0.05, na.rm=TRUE)

plotMA(res) #overexpressed genes
plotMA(resLFC) #overexpressed genes after shrinking 
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

#ONLY ERV GENES 
plotMA(res2) #NO OVEREXPRESSED GENES
plotMA(resLFC2)

#PLOT COUNTS PART
#using the overexpressed genes from TEcounts ??? 
genes <- c('Mmus38.chr1.26474418.26475755.-', 'Mmus38.chr10.71819059.71822958.-', 'Mmus38.chr12.43659073.43661217.+', 'Mmus38.chr12.49126699.49129311.+', 'Mmus38.chr12.96000997.96004887.-', 'Mmus38.chr17.41688925.41690115.-', 'Mmus38.chr17.62397216.62399306.-', 'Mmus38.chr19.48965554.48966675.+', 'Mmus38.chr19.8321592.8321936.-', 'Mmus38.chr3.49359831.49360262.+', 'Mmus38.chr7.26192088.26192735.-', 'Mmus38.chr8.97622899.97624665.+')


plotCounts(dds, gene = genes[1], intgroup = 'condition')


#WHY ALL THIS HAVE 0 counts
txi.rsem$counts[genes[1], ]
txi.rsem$counts[genes[2], ]
txi.rsem$counts[genes[3], ]
txi.rsem$counts[genes[4], ]
txi.rsem$counts[genes[5], ]
txi.rsem$counts[genes[6], ]
txi.rsem$counts[genes[7], ]
txi.rsem$counts[genes[8], ]
txi.rsem$counts[genes[9], ]
txi.rsem$counts[genes[10], ]
txi.rsem$counts[genes[11], ]
txi.rsem$counts[genes[12], ]

txi.rsem$counts['ENSMUSG00000000001', ]
txi.rsem$counts['ENSMUSG00000000028', ]
txi.rsem$abundance['ENSMUSG00000000001', ]
#DNMT3A PART
dir_path_DNMT <- '/Users/idakups/Desktop/Thesis/RsemData/DNMT3A'
samples_DNMT <- list.files(dir_path_DNMT)
metadata_DNMT <- read.table(file.path(dir_path_DNMT, 'DNMT3A_metadata.txt'), header = TRUE)

sample_paths_DNMT <- file.path(dir_path_DNMT, samples[grep(".results$", samples)])
names(sample_paths_DNMT) <- paste0("mouse", 1:9)
txi.rsem_DNMT <- tximport(sample_paths_DNMT, type = "rsem", txIn = FALSE, txOut = FALSE)

#filter out zero length transcripts
txi.rsem_DNMT$abundance <-
  txi.rsem_DNMT$abundance[apply(txi.rsem_DNMT$length,
                           1,
                           function(row) all(row !=0 )),]

txi.rsem_DNMT$counts <-
  txi.rsem_DNMT$counts[apply(txi.rsem_DNMT$length,
                        1,
                        function(row) all(row !=0 )),]

txi.rsem_DNMT$length <-
  txi.rsem_DNMT$length[apply(txi.rsem_DNMT$length,
                        1,
                        function(row) all(row !=0 )),]


#only use ERV rows
txi.rsem2_DNMT <-txi.rsem_DNMT

txi.rsem2_DNMT$abundance <- txi.rsem2_DNMT$abundance[grep("^Mmus38", rownames(txi.rsem2_DNMT$abundance)), ]
txi.rsem2_DNMT$counts <- txi.rsem2_DNMT$counts[grep("^Mmus38", rownames(txi.rsem2_DNMT$counts)), ]
txi.rsem2_DNMT$length <- txi.rsem2_DNMT$length[grep("^Mmus38", rownames(txi.rsem2_DNMT$length)), ]

print(txi.rsem_DNMT$counts)
print(dim(txi.rsem2_DNMT$counts))

# create an object storing the read counts and intermediate estimated quantities
dds_DNMT <- DESeqDataSetFromTximport(txi.rsem_DNMT, colData = metadata_DNMT, design = ~ condition  )
#dds$condition <- factor(dds$condition, levels = c("W", "M"))

dds2_DNMT <- DESeqDataSetFromTximport(txi.rsem2_DNMT, colData = metadata_DNMT, design = ~ condition  )
#dds2$condition <- factor(dds$condition, levels = c("W", "M"))

#pre-filtering - should we apply?
#smallestGroupSize <- 3 #we want to keep only if at leas 3 samples have expression greater than 10
#keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
#dds <- dds[keep,]

#STANDARD DIFFERENTIAL ANALYSIS
dds_DNMT <- DESeq(dds_DNMT)
res_DNMT <- results(dds_DNMT)
res_DNMT

dds2_DNMT <- DESeq(dds2_DNMT)
res2_DNMT <- results(dds2_DNMT)
res2_DNMT

#filter out rows with no expression
res_DNMT <- na.omit(res_DNMT)
res_DNMT

res2_DNMT <- na.omit(res2_DNMT)
res2_DNMT

#logfold change shrinkage
#here we have three condidiotns
resultsNames(dds_DNMT)
# Shrink log2 fold changes for the first coefficient
resLFC_B_vs_A <- lfcShrink(dds_DNMT, coef="condition_B_vs_A", type="apeglm")

# Shrink log2 fold changes for the second coefficient
resLFC_C_vs_A <- lfcShrink(dds_DNMT, coef="condition_C_vs_A", type="apeglm")

resultsNames(dds2_DNMT)
# Shrink log2 fold changes for the first coefficient
resLFC2_B_vs_A <- lfcShrink(dds2_DNMT, coef="condition_B_vs_A", type="apeglm")

# Shrink log2 fold changes for the second coefficient
resLFC2_C_vs_A <- lfcShrink(dds2_DNMT, coef="condition_C_vs_A", type="apeglm")




#How many adjusted p-values were lower than 0.1
sum(res_DNMT$padj < 0.1, na.rm=TRUE)
sum(res_DNMT$padj < 0.05, na.rm=TRUE)

sum(res2_DNMT$padj < 0.1, na.rm=TRUE) #5081 total
sum(res2_DNMT$padj < 0.05, na.rm=TRUE)

plotMA(res_DNMT)
plotMA(resLFC_B_vs_A)
plotMA(resLFC_C_vs_A)

plotMA(res2_DNMT)

