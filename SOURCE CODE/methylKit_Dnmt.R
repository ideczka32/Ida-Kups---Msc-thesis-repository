library(methylKit)
library(tidyverse)
library(genomation)

dir_path <- "/Users/idakups/Desktop/Thesis/methylation/DNMT3A"
dir_path
samples <- list.files(dir_path)
samples
#sample_paths <- file.path(dir_path, samples[grep("_cov.txt$", samples)])
#sample_paths <- file.path(dir_path, samples[grep("_cyt.txt$", samples)])
sample_paths <- file.path(dir_path, samples[grep(".cov$", samples)])
sample_paths
sample_paths_list <- as.list(sample_paths)
sample_paths_list

obj <- methRead(location = sample_paths_list, 
                sample.id = list( 'mut1', 'mut2','cont1', 'cont2', 'cont3'),
                assembly = 'mm10', 
                treatment = c(1, 1, 0, 0, 0),
                pipeline = "bismarkCoverage")

obj <- methRead(location = sample_paths_list, 
                sample.id = list( 'mut1', 'mut2','cont1', 'cont2', 'cont3'),
                assembly = 'mm10', 
                treatment = c(1, 1, 0, 0, 0),
                pipeline = "bismarkCytosineReport")


obj

getMethylationStats(obj[[1]],plot=FALSE,both.strands=FALSE) #sample 1
getMethylationStats(obj[[2]],plot=FALSE,both.strands=FALSE) #sample 2
getMethylationStats(obj[[3]],plot=FALSE,both.strands=FALSE) #sample 3
getMethylationStats(obj[[4]],plot=FALSE,both.strands=FALSE) #sample 4
getMethylationStats(obj[[5]],plot=FALSE,both.strands=FALSE) #sample 5

getMethylationStats(obj[[1]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj[[2]],plot=TRUE,both.strands=FALSE) #this one looks pretty bad (maybe the data qual is bad as well? i.e. some contamination as the expression is also strange in this sample)
getMethylationStats(obj[[3]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj[[4]],plot=TRUE,both.strands=FALSE) #this one looks the best
getMethylationStats(obj[[5]],plot=TRUE,both.strands=FALSE)

#All look great
getCoverageStats(obj[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[2]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[3]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[4]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[5]],plot=TRUE,both.strands=FALSE)

obj
meth <- methylKit::unite(obj, destrand=TRUE, min.per.group = 1L)
meth <- methylKit::unite(obj, min.per.group = 1L)
meth
#getCorrelation(meth,plot=TRUE)

clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

PCASamples(meth, adj.lim = c(0.5,0.1))


myDiff <- calculateDiffMeth(meth)
myDiff <- myDiff[(myDiff$chr %in% c('X', 'Y', as.character(1:22))), ]
myDiff


difference_threshold <- 25
qval_threshold <- 0.1

myDiff.hyper <- getMethylDiff(myDiff, difference=difference_threshold, qvalue = qval_threshold, type='hyper' )
myDiff.hyper #3898

myDiff.hypo <- getMethylDiff(myDiff, difference=difference_threshold, qvalue = qval_threshold, type='hypo' )
myDiff.hypo #3489


diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=qval_threshold , meth.cutoff=25)

library(genomation)
# galaxy tools converter for converting the gtf to bed 
#####USE THE MOUSE REFERENCE! 
###MAY NEED TO MAP TO NORMAL IDENTIFIERS
gene.obj=readTranscriptFeatures('/Users/idakups/Desktop/Thesis/Reference BED/mouse_full.bed',remove.unusual=FALSE) 



diffAnn <- annotateWithGeneParts(as(myDiff,"GRanges"),gene.obj)
diffAnn

diffAnn <- annotateWithGeneParts(as(myDiff.hyper,"GRanges"),gene.obj)
diffAnn

diffAnn <- annotateWithGeneParts(as(myDiff.hypo,"GRanges"),gene.obj)
diffAnn


getAssociationWithTSS <- getAssociationWithTSS(diffAnn)
getAssociationWithTSS 

association_file_path = '/Users/idakups/Desktop/Thesis/RESULTS/DNMT3A_results/associationWithTSS.txt'
write.table(getAssociationWithTSS, file=association_file_path, sep="\t", quote=FALSE, row.names=TRUE)

subset_association_with_tss_ERV <- subset(getAssociationWithTSS, grepl("^Mmus38", feature.name))
subset_association_with_tss_ERV
n_row_ERV <- nrow(subset_association_with_tss_ERV)
n_row_ERV #66232

subset_association_with_tss_gene <- subset(getAssociationWithTSS, grepl("^ENS", feature.name))
subset_association_with_tss_gene
n_row_gene <- nrow(subset_association_with_tss_gene)
n_row_gene #908444


#ONLY COUNT METHYLATION AROUND PROMOTER
subset_association_with_tss_ERV <- subset(subset_association_with_tss_ERV, abs(dist.to.feature) < 2000)
n_row_ERV <- nrow(subset_association_with_tss_ERV)
n_row_ERV #16264

subset_association_with_tss_gene <- subset(subset_association_with_tss_gene, abs(dist.to.feature) < 2000)
n_row_gene <- nrow(subset_association_with_tss_gene)
n_row_gene #677905

unique_feature_names_erv <- unique(subset_association_with_tss_ERV$feature.name)
num_unique_feature_names_erv <- length(unique_feature_names_erv)
num_unique_feature_names_erv #5426
x<- table(subset_association_with_tss_ERV$feature.name )
count_greater_than_5 <- sum(x > 30)
count_greater_than_5
names_greater_than_50 <- names(x[x > 30])
names_greater_than_50


unique_feature_names_gene <- unique(subset_association_with_tss_gene$feature.name)
num_unique_feature_names_gene <- length(unique_feature_names_gene)
num_unique_feature_names_gene #54453
x<- table(subset_association_with_tss_gene$feature.name )
count_greater_than_5 <- sum(x > 30)
count_greater_than_5

