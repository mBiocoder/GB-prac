#BiocManager::install("DEXSeq")
#BiocManager::install("GenomicRanges")
#install.packages("ROCR")
library(DEXSeq)
library(GenomicRanges)
library(ROCR)
library(tidyr)


############### generate DEXSeq results ########################################
## list the htseq-count files
countFiles = list.files("extdata/", full=T)
names(countFiles) <- gsub("_dexseq.txt", "", basename(countFiles))

## prepare the sample annotation
group <- rep(1,length(countFiles))
# sample names ending in 6..10 are in group 2
group[grepl("[06-9]$", names(countFiles))] <- 2
sampleTable = data.frame(condition=factor(group))
rownames(sampleTable) = names(countFiles)

## load the data
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile = "annotation_dexseq_b37.gff" )

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



################## analysis ####################################################
# data from our analysis (generated before)
our_result <- read.delim("diff_analysis_results.tsv")
our_result <- separate(our_result, exon, c("exon_start", "exon_end"), "-")
our_result$exon_start <- as.numeric(our_result$exon_start)
our_result$exon_end <- as.numeric(our_result$exon_end)-1
our_result <- our_result[, c(1,2,3,10,11)]
exon <- paste(our_result$gene, paste(our_result$exon_start,our_result$exon_end,sep="-"), sep=":")
our_result <- cbind(our_result, exon)

# data from the DEXSeq analysis
dexSeq_result <- read.delim("dexSeq_results.tsv")
colnames(dexSeq_result)[9] <- "exon_start"
colnames(dexSeq_result)[10] <- "exon_end"
colnames(dexSeq_result)[1] <- "gene"
dexSeq_result <- separate_rows(dexSeq_result, 1, sep = "\\+")
exon <- paste(dexSeq_result$gene, paste(dexSeq_result$exon_start,dexSeq_result$exon_end,sep="-"), sep=":")
dexSeq_result <- cbind(dexSeq_result, exon)
# filter out rows with exonBaseMean = 0.000
dexSeq_result <- dexSeq_result[dexSeq_result$exonBaseMean != 0.0, ]
dexSeq_result <- dexSeq_result[, c(1,6,7,9,10,24)]


# load the ground truth data
load("differential_exons.RData")
ground_truth <- data.frame(differential.skipped)
ground_truth <- ground_truth[,c(10, 2, 3, 7)]
ground_truth <- ground_truth[ground_truth$type == "exon",]
ground_truth <- ground_truth[, c(1,2,3)]
colnames(ground_truth) <- c("gene", "exon_start", "exon_end")
exon <- paste(ground_truth$gene, paste(ground_truth$exon_start,ground_truth$exon_end,sep="-"), sep=":")
ground_truth <- cbind(ground_truth, exon)


## ROC-plots ###
our_labels <- as.numeric(our_result$exon %in% ground_truth$exon)
pred <- prediction(our_result$padj, our_labels, label.ordering = c(1,0))
perf <- performance(pred, "tpr", "fpr")
plot(perf, col="red", main="ROC curve of differentially spliced exons")


# choose only the exons from the DEXSeq results that are also in the LRS result
dexSeq_result <- merge(our_result[, c(1,2,3,6)], dexSeq_result, by=c("gene", "exon_start", "exon_end", "exon"))
# remove (padj) NA lines
dexSeq_result <- dexSeq_result[!is.na(dexSeq_result$padj),]

dex_labels <- as.numeric(dexSeq_result$exon %in% ground_truth$exon)
dex_pred <- prediction(dexSeq_result$padj, dex_labels, label.ordering = c(1,0))
dex_perf <- performance(dex_pred, "tpr", "fpr")
plot(dex_perf, col="blue", add=TRUE)

### AUC #######
lrs.auc.perf = performance(pred, measure = "auc")
lrs.auc.perf@y.values
dex.auc.perf = performance(dex_pred, measure = "auc")
dex.auc.perf@y.values

lrsLab <- paste0("LRS (AUC: ", paste0(round(lrs.auc.perf@y.values[[1]], 4), ")"))
dexLab <- paste0("DEXSeq (AUC: ", paste0(round(dex.auc.perf@y.values[[1]], 4), ")"))

legend("bottomright", legend = c(lrsLab, dexLab), col = c("red", "blue"), pch = 19, bty = "n")
abline(a=0, b= 1)


## compute sensitivity and FDR #################################################
# sensitivity = TP / (TP + FN) fraction of correctly classified exons out of all differential exons
# FDR = FP / (FP + TP) fraction of not differential exons out of all exons classified as differential

# TP: prediction is positive (1) and truth positive (1)
# FN: prediction is negative (0) and truth positive (1)
# FP: prediction is positive (1) and truth negative (0)

# use binary vectors with 0/FALSE = not differentially spliced/padj >= cutoff
# 1/TRUE = differentially spliced/padj < cutoff

sensitivity <- function(cutoff, padj, labels){
  prediction <- as.numeric(padj < cutoff)
  TP <- sum(labels & prediction)
  FN <- sum(labels & !prediction)
  print("TP")
  print(TP)
  print("FN")
  print(FN)
  sens <- TP / (TP+FN)
  return(sens)
}

FDR <- function(cutoff, padj, labels){
  prediction <- as.numeric(padj < cutoff)
  FP <- sum(!labels & prediction)
  TP <- sum(labels & prediction)
  print("FP")
  print(FP)
  print("TP")
  print(TP)
  fdr <- FP / (FP+TP)
  return(fdr)
}


sensitivity(0.05, our_result$padj, our_labels)
FDR(0.05, our_result$padj, our_labels)

sensitivity(0.05, dexSeq_result$padj, dex_labels)
FDR(0.05, dexSeq_result$padj, dex_labels)