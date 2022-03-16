library(dplyr)
library(DEXSeq)

source("load_SubreadOutput.R")
samp <- data.frame(row.names = c("data_all_samples_STAR/WT_0h_rep_1_star.bam", "data_all_samples_STAR/WT_0h_rep_2_star.bam",
"data_all_samples_STAR/WT_0h_rep_3_star.bam", "data_all_samples_STAR/WT_0h_rep_4_star.bam", "data_all_samples_STAR/WT_1h_rep_1_star.bam",
"data_all_samples_STAR/WT_1h_rep_2_star.bam", "data_all_samples_STAR/WT_1h_rep_3_star.bam", "data_all_samples_STAR/WT_1h_rep_4_star.bam",
"data_all_samples_STAR/WT_2h_rep_1_star.bam", "data_all_samples_STAR/WT_2h_rep_2_star.bam", "data_all_samples_STAR/WT_2h_rep_3_star.bam",
"data_all_samples_STAR/WT_2h_rep_4_star.bam", "data_all_samples_STAR/WT_4h_rep_1_star.bam", "data_all_samples_STAR/WT_4h_rep_2_star.bam",
"data_all_samples_STAR/WT_4h_rep_3_star.bam", "data_all_samples_STAR/WT_4h_rep_4_star.bam"), 
                   condition = rep(c("WT_0h","WT_1h", "WT_2h", "WT_4h"),each=4))


dxd <- DEXSeqDataSetFromFeatureCounts("Mus_musculus_38_75_featureCount.out",
                                         flattenedfile = "Mus_musculus.GRCm38.75_flat.gtf",sampleData = samp)

## estimate the size factors
dxd = estimateSizeFactors( dxd ) 
## fit function to compute the variance in dependence on the mean
dxd = estimateDispersions( dxd ) 
## check the fit of the variance (dispersion) model
plotDispEsts( dxd ) 

## run the test for differential exon usage
dxd = testForDEU( dxd ) 
## obtain the results
dxr1 = DEXSeqResults( dxd )
# save the result so you don't always have to repeat DEXSeq
write.table(dxr1, file="dexSeq_results.tsv", sep='\t', quote=FALSE, row.names = FALSE)

#read dexseq results
dxr1 <- read.delim("./Resulting_tables/Set_of_alternatively_spliced_genes.tsv")

#make copy and remove all NAs
dxr_copy <- dxr1
dxr_copy <- na.omit(dxr_copy)

#explore results
mcols(dxr_copy)$description
table ( dxr_copy$padj < 0.01 )

#visualization
dxr2 <- dxr_copy[, dxr_copy$padj < 0.01]

#Vegfa gene: ENSMUSG00000023951
plotDEXSeq( dxr1, "ENSMUSG00000023951", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "ENSMUSG00000023951", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "ENSMUSG00000023951", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


#Hmgcr
plotDEXSeq( dxr1, "ENSMUSG00000021670", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr1, "ENSMUSG00000021670", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr1, "ENSMUSG00000021670", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#Apobec1
plotDEXSeq( dxr1, "ENSMUSG00000040613", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr1, "ENSMUSG00000040613", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#Apobr !! -> Paper
plotDEXSeq( dxr1, "ENSMUSG00000042759", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr1, "ENSMUSG00000042759", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#lox1

plotDEXSeq( dxr1, "ENSMUSG00000024529", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr1, "ENSMUSG00000040613", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#foxp4
plotDEXSeq( dxr1, "ENSMUSG00000023991", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr1, "ENSMUSG00000023991", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#foxp1
plotDEXSeq( dxr1, "ENSMUSG00000030067+ENSMUSG00000030068+ENSMUSG00000093661", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr1, "ENSMUSG00000030067", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

