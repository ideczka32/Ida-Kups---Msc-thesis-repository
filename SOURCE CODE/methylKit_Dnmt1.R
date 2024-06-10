library(methylKit)
library(tidyverse)
library(genomation)


dir_path <- "/Users/idakups/Desktop/Thesis/methylation/DNMT1"
samples <- list.files(dir_path)
sample_paths <- file.path(dir_path, samples[grep("cov$", samples)])
sample_paths_list <- as.list(sample_paths)

# Now use file_list to create the object
obj <- methRead(location = sample_paths_list, 
                sample.id = list('cont1', 'cont2', 'mut1', 'mut2'),
                assembly = 'hg38', 
                treatment = c(0, 0, 1, 1),
                pipeline = "bismarkCoverage")

getMethylationStats(obj[[1]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj[[2]],plot=TRUE,both.strands=FALSE) #this one looks pretty bad (maybe the data qual is bad as well? i.e. some contamination as the expression is also strange in this sample)
getMethylationStats(obj[[3]],plot=TRUE,both.strands=FALSE)
getMethylationStats(obj[[4]],plot=TRUE,both.strands=FALSE)


#All look great
getCoverageStats(obj[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[2]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[3]],plot=TRUE,both.strands=FALSE)
getCoverageStats(obj[[4]],plot=TRUE,both.strands=FALSE) 

#should i sort the files first? 
#could use min.per.group to include the reads covered only in one member of the group 
meth <- methylKit::unite(obj) #, destrand=TRUE
meth


clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

PCASamples(meth)

#CONSIDER TILING WINDOW ANALYSIS

#DIFFERENTIAL METHYLATION 
myDiff <- calculateDiffMeth(meth)

difference_threshold <- 25
qval_threshold <- 0.1

myDiff.hyper <- getMethylDiff(myDiff, difference=difference_threshold, qvalue = qval_threshold, type='hyper' )
myDiff.hyper #12427

myDiff.hypo <- getMethylDiff(myDiff, difference=difference_threshold, qvalue = qval_threshold, type='hypo' )
myDiff.hypo #2136565

myDiff

# Subset the data to filter out rows where 'chr' is not 'X', 'Y', or any number from 1 to 22
filtered_data <- myDiff[!(myDiff$chr %in% c('X', 'Y', as.character(1:22))), ]

# Get unique values in 'chr' column after filtering
unique_chr_values <- unique(filtered_data$chr)
unique_chr_values


diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=qval_threshold , meth.cutoff=meth_cutoff) #CYTOSINE REPORT
diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=qval_threshold , meth.cutoff=meth_cutoff, exclude = unique_chr_values)


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

association_file_path = '/Users/idakups/Desktop/Thesis/RESULTS/DNMT1_results/associationWithTSS.txt'
write.table(getAssociationWithTSS, file=association_file_path, sep="\t", quote=FALSE, row.names=TRUE)

subset_association_with_tss_ERV <- subset(getAssociationWithTSS, grepl("^Hsap38", feature.name))
subset_association_with_tss_ERV
n_row_ERV <- nrow(subset_association_with_tss_ERV)
n_row_ERV #259441

subset_association_with_tss_gene <- subset(getAssociationWithTSS, grepl("^ENS", feature.name))
subset_association_with_tss_gene
n_row_gene <- nrow(subset_association_with_tss_gene)
n_row_gene #581053


#ONLY COUNT METHYLATION AROUND PROMOTER
subset_association_with_tss_ERV <- subset(subset_association_with_tss_ERV, abs(dist.to.feature) < 2000)
n_row_ERV <- nrow(subset_association_with_tss_ERV)
n_row_ERV #27431

subset_association_with_tss_gene <- subset(subset_association_with_tss_gene, abs(dist.to.feature) < 2000)
n_row_gene <- nrow(subset_association_with_tss_gene)
n_row_gene #100701

unique_feature_names_erv <- unique(subset_association_with_tss_ERV$feature.name)
num_unique_feature_names_erv <- length(unique_feature_names_erv)
num_unique_feature_names_erv 
x<- table(subset_association_with_tss_ERV$feature.name )
count_greater_than_50 <- sum(x > 50)
count_greater_than_50 #89 

names_greater_than_50 <- names(x[x > 50])
names_greater_than_50

write.table(names_greater_than_50, file = "/Users/idakups/Desktop/Thesis/RESULTS/DNMT1_results/names_greater_than_30.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


unique_feature_names_gene <- unique(subset_association_with_tss_gene$feature.name)
num_unique_feature_names_gene <- length(unique_feature_names_gene)
num_unique_feature_names_gene 
x<- table(subset_association_with_tss_gene$feature.name )
count_greater_than_50 <- sum(x > 50)
count_greater_than_50 #348
