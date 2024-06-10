#TET2
library(methylKit)
library(tidyverse)
library(genomation)

dir_path <- "/Users/idakups/Desktop/Thesis/methylation/TET2"
samples <- list.files(dir_path)
sample_paths <- file.path(dir_path, samples[grep("_cyt.txt$", samples)])
sample_cov_paths <- file.path(dir_path, samples[grep("_cov.txt$", samples)])
sample_paths_list <- as.list(sample_paths)
sample_cov_paths_list <- as.list(sample_cov_paths)

sample_cov_paths_list_sinmut2 <- sample_cov_paths_list[-(length(sample_cov_paths_list) - 1)]


# Now use file_list to create the object
obj <- methRead(location = sample_paths_list, 
                sample.id = list('cont1', 'cont2', 'mut1', 'mut2', 'mut3'),
                assembly = 'mm10', 
                treatment = c(0, 0, 1, 1, 1),
                pipeline = "bismarkCytosineReport")


obj_cov <- methRead(location = sample_cov_paths_list, 
                sample.id = list('cont1', 'cont2', 'mut1', 'mut2', 'mut3'),
                assembly = 'mm10', 
                treatment = c(0, 0, 1, 1, 1),
                pipeline = "bismarkCoverage")

obj_cov_sinmut2 <- methRead(location = sample_cov_paths_list_sinmut2, 
                            sample.id = list('cont1', 'cont2', 'mut1', 'mut3'),
                            assembly = 'mm10', 
                            treatment = c(0, 0, 1, 1),
                            pipeline = "bismarkCoverage")


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
getCoverageStats(obj[[4]],plot=TRUE,both.strands=FALSE) #REALLY LOW READ COVERAGE
getCoverageStats(obj[[5]],plot=TRUE,both.strands=FALSE)

getMethylationStats(obj_cov[[1]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj_cov[[2]],plot=TRUE,both.strands=FALSE) #this one looks pretty bad (maybe the data qual is bad as well? i.e. some contamination as the expression is also strange in this sample)
getMethylationStats(obj_cov[[3]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj_cov[[4]],plot=TRUE,both.strands=FALSE) #this one looks the best
getMethylationStats(obj_cov[[5]],plot=TRUE,both.strands=FALSE)

#All look great
getCoverageStats(obj_cov[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj_cov[[2]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj_cov[[3]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj_cov[[4]],plot=TRUE,both.strands=FALSE) #REALLY LOW READ COVERAGE
getCoverageStats(obj_cov[[5]],plot=TRUE,both.strands=FALSE)

# creates a methylBase object, 
# where only CpGs covered with at least 1 sample per group will be returned

# there were two groups defined by the treatment vector, 
# given during the creation of myobj: treatment=c(1,1,0,0)
#meth.min=unite(obj,min.per.group=1L)


#get the cytosines covered in all samples min.per.group
meth <- methylKit::unite(obj, destrand=TRUE)
meth <- methylKit::unite(obj, min.per.group = 1L) #this is for coverage object
meth

ERV_df <- meth[str_starts(rownames(res), 'Mmus38'), ]

#save to file 
file_path = '/Users/idakups/Desktop/Thesis/RsemData/TET2/TET2_methylkit_erv_full_results.tsv'
write.table(ERV_df, file = file_path, sep = "\t", row.names = TRUE)



clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

PCASamples(meth)

#COV sinmut2
meth_cov_sinmut2 <- methylKit::unite(obj_cov_sinmut2)
meth_cov_sinmut2

clusterSamples(meth_cov_sinmut2, dist="correlation", method="ward", plot=TRUE)
clusterSamples(meth_cov, dist="correlation", method="ward", plot=TRUE)

PCASamples(meth_cov)

#cov 
meth_cov <- methylKit::unite(obj_cov)
meth_cov

meth_covx <- methylKit::unite(obj_cov, min.per.group = 1L)
meth_covx



clusterSamples(meth_cov, dist="correlation", method="ward", plot=TRUE)

PCASamples(meth_cov)



#DIFFERENTIAL METHYLATION
myDiff <- calculateDiffMeth(meth)
diff.threshold <- 20

myDiff_cov <- calculateDiffMeth(meth_cov)
myDiff_covx <- calculateDiffMeth(meth_covx)

myDiff_sinmut2 <- calculateDiffMeth(meth_cov_sinmut2)

#GET HYPERMETHYLATED BASES WITH DIFF > 25%
myDiff25p.hyper=getMethylDiff(myDiff,difference=diff.threshold,qvalue=0.,type="hyper")
myDiff25p.hyper_cov=getMethylDiff(myDiff_cov,difference=diff.threshold,qvalue=0.1,type="hyper")

#GET HYPOMETHYLATED BASES WITH DIFF > 25%
myDiff25p.hypo=getMethylDiff(myDiff,difference=diff.threshold,qvalue=0.1,type="hypo")
myDiff25p.hypo_cov=getMethylDiff(myDiff_cov,difference=diff.threshold,qvalue=0.1,type="hypo")

#GET ALL DIFFERENTIALLY METHYLATED BASES
myDiff25p=getMethylDiff(myDiff,difference=diff.threshold,qvalue=0.1)


myDiff25p_cov=getMethylDiff(myDiff_cov,difference=diff.threshold,qvalue=0.1)
myDiff25p_covx=getMethylDiff(myDiff_covx,difference=diff.threshold,qvalue=0.1)


myDiff25p
myDiff25p_cov
myDiff25p_covx

meth_cutoff <- 25 
diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.2, meth.cutoff=meth_cutoff) #CYTOSINE REPORT
diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.2, meth.cutoff=meth_cutoff)

diffMethPerChr(myDiff_cov,plot=FALSE,qvalue.cutoff=0.2, meth.cutoff=meth_cutoff)
diffMethPerChr(myDiff_cov,plot=TRUE,qvalue.cutoff=0.2, meth.cutoff=meth_cutoff) #COVERAGE  
diffMethPerChr(myDiff_covx,plot=TRUE,qvalue.cutoff=0.2, meth.cutoff=meth_cutoff) 

diffMethPerChr(myDiff_sinmut2,plot=FALSE,qvalue.cutoff=0.2, meth.cutoff=meth_cutoff)
diffMethPerChr(myDiff_sinmut2,plot=TRUE,qvalue.cutoff=0.2, meth.cutoff=meth_cutoff) #COVERAGE  
#ERV's differentialy expressed on chrom 
#[15,17,4,8]
#ERV's expressed on this chrom
#['1' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '2' '3' '4' '5' '6' '7' '8' '9' 'X' 'Y']
#ERV's sum expr > 50
#['1' '11' '12' '15' '16' '17' '18' '19' '2' '3' '4' '5' '6' '7' '8' '9''X']
#ERV's sum expr > 100
#['1' '11' '12' '16' '17' '18' '19' '3' '4' '5' '6' '7' '8' 'X']

#BETTER AVOID USING STATS
chr15 <- myDiff25p[myDiff25p$chr == 15]



chr17 <- myDiff25p[myDiff25p$chr == 17]
chr4 <- myDiff25p[myDiff25p$chr == 4]
chr8 <- myDiff25p[myDiff25p$chr == 8]

by_chrom <- split(meth, meth$chr)
by_chrom

write.table(meth, file='/Users/idakups/Desktop/Thesis/methylation/TET2/meth_TET2.tsv', quote=FALSE, sep='\t')
write.table(meth_covx, file='/Users/idakups/Desktop/Thesis/methylation/TET2/meth_cov_TET2.tsv', quote=FALSE, sep='\t')


#ANNOTATING DIFFERENTIALLY METHYLATED GENES/REGIONS GENOMATION
library(genomation)

gene.obj=readTranscriptFeatures('/Users/idakups/Desktop/Thesis/Reference BED/mouse_full.bed',remove.unusual=FALSE)
diffAnn <- annotateWithGeneParts(as(myDiff25p_covx,"GRanges"),gene.obj)
diffAnn2 <- annotateWithGeneParts(as(myDiff25p_covx,"GRanges"),gene.obj,intersect.chr=TRUE )
diffAnnn
diffAnn2

gene.obj=readTranscriptFeatures('/Users/idakups/Desktop/Thesis/Reference BED/mouse_full.bed',remove.unusual=FALSE)
diffAnn <- annotateWithGeneParts(as(myDiff,"GRanges"),gene.obj)



getAssociationWithTSS <- getAssociationWithTSS(diffAnn)
getAssociationWithTSS #fom my 160 differentially methylated regions none has the feature name starting with Mmus

subset_association_with_tss_ERV <- subset(getAssociationWithTSS, grepl("^Mmus38", feature.name))
subset_association_with_tss_ERV
n_row_ERV <- nrow(subset_association_with_tss_ERV)
n_row_ERV #68

subset_association_with_tss_gene <- subset(getAssociationWithTSS, grepl("^ENS", feature.name))
subset_association_with_tss_gene
n_row_gene <- nrow(subset_association_with_tss_gene)
n_row_gene #150


#ONLY COUNT METHYLATION AROUND PROMOTER
subset_association_with_tss_ERV <- subset(subset_association_with_tss_ERV, abs(dist.to.feature) < 2000)
n_row_ERV <- nrow(subset_association_with_tss_ERV)
n_row_ERV #14841

subset_association_with_tss_gene <- subset(subset_association_with_tss_gene, abs(dist.to.feature) < 2000)
n_row_gene <- nrow(subset_association_with_tss_gene)
n_row_gene #569249

unique_feature_names_erv <- unique(subset_association_with_tss_ERV$feature.name)
num_unique_feature_names_erv <- length(unique_feature_names_erv)
num_unique_feature_names_erv #5426
x<- table(subset_association_with_tss_ERV$feature.name )
count_greater_than_5 <- sum(x > 2)
count_greater_than_5


unique_feature_names_gene <- unique(subset_association_with_tss_gene$feature.name)
num_unique_feature_names_gene <- length(unique_feature_names_gene)
num_unique_feature_names_gene #54453
x<- table(subset_association_with_tss_gene$feature.name )
count_greater_than_5 <- sum(x > 30)
count_greater_than_5








