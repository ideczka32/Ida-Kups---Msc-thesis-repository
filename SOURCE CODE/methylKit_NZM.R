library(methylKit)
library(tidyverse)
dir_path <- "/Users/idakups/Desktop/Thesis/methylation/NZM"
samples <- list.files(dir_path)


sample_paths <- file.path(dir_path, samples[grep("_cyt.txt$", samples)])
sample_paths <- file.path(dir_path, samples[grep("cov$", samples)])
sample_paths_list <- as.list(sample_paths)

obj <- methRead(location = sample_paths_list, 
                sample.id = list('inv1', 'inv2','inv3', 'n-inv1', 'n-inv2', 'inv4', 'inv5', 'n-inv3', 'n-inv4'),
                assembly = 'hg38', 
                treatment = c(1, 1, 1, 0, 0,1,1,0,0),
                pipeline = "bismarkCytosineReport")

obj <- methRead(location = sample_paths_list, 
                sample.id = list('inv1', 'inv2','inv3', 'n-inv1', 'n-inv2', 'inv4', 'inv5', 'n-inv3', 'n-inv4'),
                assembly = 'hg38', 
                treatment = c(1, 1, 1, 0, 0,1,1,0,0),
                context = 'CpG',
                pipeline = "bismarkCoverage",
                sep='\t')


obj
getMethylationStats(obj[[1]],plot=FALSE,both.strands=FALSE) #sample 1
getMethylationStats(obj[[2]],plot=FALSE,both.strands=FALSE) #sample 2
getMethylationStats(obj[[3]],plot=FALSE,both.strands=FALSE) #sample 3
getMethylationStats(obj[[4]],plot=FALSE,both.strands=FALSE) #sample 4
getMethylationStats(obj[[5]],plot=FALSE,both.strands=FALSE) #sample 5

#all plots look quite okay 
getMethylationStats(obj[[1]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj[[2]],plot=TRUE,both.strands=FALSE) 
getMethylationStats(obj[[3]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj[[4]],plot=TRUE,both.strands=FALSE) 
getMethylationStats(obj[[5]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj[[6]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj[[7]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj[[8]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj[[9]],plot=TRUE,both.strands=FALSE)

#All look great
getCoverageStats(obj[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[2]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[3]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[4]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[5]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[6]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[7]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[8]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[9]],plot=TRUE,both.strands=FALSE)


meth <- methylKit::unite(obj, min.per.group = 1L)
meth


#getCorrelation(meth,plot=TRUE)

clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

PCASamples(meth, adj.lim = c(0.5,0.1))

myDiff <- calculateDiffMeth(meth)
filtered_data <- myDiff[(myDiff$chr %in% c('X', 'Y', as.character(1:22))), ]
filtered_data


diffMethPerChr(filtered_data,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=20)

myDiff.hyper <- getMethylDiff(filtered_data, difference=difference_threshold, qvalue = qval_threshold, type='hyper' )
myDiff.hyper #188907

myDiff.hypo <- getMethylDiff(filtered_data, difference=difference_threshold, qvalue = qval_threshold, type='hypo' )
myDiff.hypo #231788


#ANNOTATING DIFFERENTIALLY METHYLATED GENES/REGIONS GENOMATION
library(genomation)
#using gtf2bed to obtain the bed file 
########THIS ONE IS NOT WORKING
# galaxy tools converter for converting the gtf to bed 
#gene.obj=readTranscriptFeatures('/Users/idakups/Desktop/Thesis/Reference BED/h.bed',remove.unusual=FALSE) 


gene.obj.hsap=readTranscriptFeatures('/Users/idakups/Desktop/Thesis/Reference BED/human_full.bed',remove.unusual=FALSE) #Doesnt work probably because of the erv naming format ...


diffAnn <- annotateWithGeneParts(as(filtered_data,"GRanges"),gene.obj.hsap)


getAssociationWithTSS <- getAssociationWithTSS(diffAnn)
getAssociationWithTSS #regions ng with Mmus

subset_association_with_tss_ERV <- subset(getAssociationWithTSS, grepl("^Hsap38", feature.name))
subset_association_with_tss_ERV
n_row_ERV <- nrow(subset_association_with_tss_ERV)
n_row_ERV #229914

subset_association_with_tss_gene <- subset(getAssociationWithTSS, grepl("^ENS", feature.name))
subset_association_with_tss_gene
n_row_gene <- nrow(subset_association_with_tss_gene)
n_row_gene # 3728474

subset_association_with_tss_ERV <- subset(subset_association_with_tss_ERV, abs(dist.to.feature) < 2000)
n_row_ERV <- nrow(subset_association_with_tss_ERV)
n_row_ERV #21043

subset_association_with_tss_gene <- subset(subset_association_with_tss_gene, abs(dist.to.feature) < 2000)
n_row_gene <- nrow(subset_association_with_tss_gene)
n_row_gene #2339711

unique_feature_names_erv <- unique(subset_association_with_tss_ERV$feature.name)
num_unique_feature_names_erv <- length(unique_feature_names_erv)
num_unique_feature_names_erv 
x<- table(subset_association_with_tss_ERV$feature.name )
count_greater_than_50 <- sum(x > 50)
count_greater_than_50 # thr 30 54 and thr 50 8 

names_greater_than_50 <- names(x[x > 50])
names_greater_than_50




