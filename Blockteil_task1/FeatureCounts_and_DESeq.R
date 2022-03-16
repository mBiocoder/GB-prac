library(Rsubread)
library(Rsamtools)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(pander)
library(KEGGGraph)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(GO.db)
library(topGO)
library(dplyr)
library(pathview)
library(org.Mm.eg.db)
library(limma)
library(edgeR)
library(VennDiagram)
library(ImpulseDE2)
library(GenomicFeatures)
library(HTSFilter)
library(SPIA)
library(DESeq2)
library(DOSE)
library(enrichplot)


##################################### Run featureCounts #######################################

#### STAR ####
#Generate single count file containing all samples
file.list <- list.files("./data_all_samples_STAR")

all_samples.output.bamFile <- c("./data_all_samples_STAR/WT_0h_rep_1_star.bam", "./data_all_samples_STAR/WT_0h_rep_2_star.bam",
                                "./data_all_samples_STAR/WT_0h_rep_3_star.bam", "./data_all_samples_STAR/WT_0h_rep_4_star.bam",
                                "./data_all_samples_STAR/WT_1h_rep_1_star.bam", "./data_all_samples_STAR/WT_1h_rep_2_star.bam",
                                "./data_all_samples_STAR/WT_1h_rep_3_star.bam", "./data_all_samples_STAR/WT_1h_rep_4_star.bam",
                                "./data_all_samples_STAR/WT_2h_rep_1_star.bam", "./data_all_samples_STAR/WT_2h_rep_2_star.bam",
                                "./data_all_samples_STAR/WT_2h_rep_3_star.bam", "./data_all_samples_STAR/WT_2h_rep_4_star.bam",
                                "./data_all_samples_STAR/WT_4h_rep_1_star.bam", "./data_all_samples_STAR/WT_4h_rep_2_star.bam",
                                "./data_all_samples_STAR/WT_4h_rep_3_star.bam", "./data_all_samples_STAR/WT_4h_rep_4_star.bam")


all_samples_counts <- featureCounts(all_samples.output.bamFile, 
                                    annot.inbuilt="mm10",
                                    isPairedEnd=TRUE)

#Write to output
RESULTS_DIR <- "./results/"
write.table(all_samples_counts$counts, file=file.path(RESULTS_DIR,"all_samples_raw_read_counts.txt"), 
            sep="\t", quote=F, append=F)



#### hisat ####
#Generate single count file containing all samples
file.list.hisat <- list.files("./all_samples_hisat_bam")

all_samples.output.bamFile.hisat <- c("./all_samples_hisat_bam/WT_0h_rep_1_hisat.bam", "./all_samples_hisat_bam/WT_0h_rep_2_hisat.bam",
                                      "./all_samples_hisat_bam/WT_0h_rep_3_hisat.bam", "./all_samples_hisat_bam/WT_0h_rep_4_hisat.bam",
                                      "./all_samples_hisat_bam/WT_1h_rep_1_hisat.bam", "./all_samples_hisat_bam/WT_1h_rep_2_hisat.bam",
                                      "./all_samples_hisat_bam/WT_1h_rep_3_hisat.bam", "./all_samples_hisat_bam/WT_1h_rep_4_hisat.bam",
                                      "./all_samples_hisat_bam/WT_2h_rep_1_hisat.bam", "./all_samples_hisat_bam/WT_2h_rep_2_hisat.bam",
                                      "./all_samples_hisat_bam/WT_2h_rep_3_hisat.bam", "./all_samples_hisat_bam/WT_2h_rep_4_hisat.bam",
                                      "./all_samples_hisat_bam/WT_4h_rep_1_hisat.bam", "./all_samples_hisat_bam/WT_4h_rep_2_hisat.bam",
                                      "./all_samples_hisat_bam/WT_4h_rep_3_hisat.bam", "./all_samples_hisat_bam/WT_4h_rep_4_hisat.bam")


all_samples_counts_hisat <- featureCounts(all_samples.output.bamFile.hisat, 
                                          annot.inbuilt="mm10",
                                          isPairedEnd=TRUE)

#Write to output
RESULTS_DIR <- "./results/hisat/"
write.table(all_samples_counts_hisat$counts, file=file.path(RESULTS_DIR,"hisat_all_samples_raw_read_counts.txt"), 
            sep="\t", quote=F, append=F)


#### contextmap ####
#Generate single count file containing all samples
file.list.contextmap <- list.files("./all_samples_contextmap_bam")

all_samples.output.bamFile.contextmap <- c("./all_samples_contextmap_bam/WT_0h_rep_1_contextmap.bam", "./all_samples_contextmap_bam/WT_0h_rep_2_contextmap.bam",
                                           "./all_samples_contextmap_bam/WT_0h_rep_3_contextmap.bam", "./all_samples_contextmap_bam/WT_0h_rep_4_contextmap.bam",
                                           "./all_samples_contextmap_bam/WT_1h_rep_1_contextmap.bam", "./all_samples_contextmap_bam/WT_1h_rep_2_contextmap.bam",
                                           "./all_samples_contextmap_bam/WT_1h_rep_3_contextmap.bam", "./all_samples_contextmap_bam/WT_1h_rep_4_contextmap.bam",
                                           "./all_samples_contextmap_bam/WT_2h_rep_1_contextmap.bam", "./all_samples_contextmap_bam/WT_2h_rep_2_contextmap.bam",
                                           "./all_samples_contextmap_bam/WT_2h_rep_3_contextmap.bam", "./all_samples_contextmap_bam/WT_2h_rep_4_contextmap.bam",
                                           "./all_samples_contextmap_bam/WT_4h_rep_1_contextmap.bam", "./all_samples_contextmap_bam/WT_4h_rep_2_contextmap.bam",
                                           "./all_samples_contextmap_bam/WT_4h_rep_3_contextmap.bam", "./all_samples_contextmap_bam/WT_4h_rep_4_contextmap.bam")


all_samples_counts_contextmap <- featureCounts(all_samples.output.bamFile.contextmap, 
                                               annot.inbuilt="mm10",
                                               isPairedEnd=TRUE)

#Write to output
RESULTS_DIR <- "./results/contextmap/"
write.table(all_samples_counts_contextmap$counts, file=file.path(RESULTS_DIR,"contextmap_all_samples_raw_read_counts.txt"), 
            sep="\t", quote=F, append=F)


#### tophat2 ####
#Generate single count file containing all samples
file.list.tophat2 <- list.files("./all_samples_tophat2_bam")

all_samples.output.bamFile.tophat2 <- c("./all_samples_tophat2_bam/WT_0h_rep_1_tophat2.bam", "./all_samples_tophat2_bam/WT_0h_rep_2_tophat2.bam",
                                        "./all_samples_tophat2_bam/WT_0h_rep_3_tophat2.bam", "./all_samples_tophat2_bam/WT_0h_rep_4_tophat2.bam",
                                        "./all_samples_tophat2_bam/WT_1h_rep_1_tophat2.bam", "./all_samples_tophat2_bam/WT_1h_rep_2_tophat2.bam",
                                        "./all_samples_tophat2_bam/WT_1h_rep_3_tophat2.bam", "./all_samples_tophat2_bam/WT_1h_rep_4_tophat2.bam",
                                        "./all_samples_tophat2_bam/WT_2h_rep_1_tophat2.bam", "./all_samples_tophat2_bam/WT_2h_rep_2_tophat2.bam",
                                        "./all_samples_tophat2_bam/WT_2h_rep_3_tophat2.bam", "./all_samples_tophat2_bam/WT_2h_rep_4_tophat2.bam",
                                        "./all_samples_tophat2_bam/WT_4h_rep_1_tophat2.bam", "./all_samples_tophat2_bam/WT_4h_rep_2_tophat2.bam",
                                        "./all_samples_tophat2_bam/WT_4h_rep_3_tophat2.bam", "./all_samples_tophat2_bam/WT_4h_rep_4_tophat2.bam")


all_samples_counts_tophat2 <- featureCounts(all_samples.output.bamFile.tophat2, 
                                            annot.inbuilt="mm10",
                                            isPairedEnd=TRUE)

#Write to output
RESULTS_DIR <- "./results/tophat2/"
write.table(all_samples_counts_tophat2$counts, file=file.path(RESULTS_DIR,"tophat2_all_samples_raw_read_counts.txt"), 
            sep="\t", quote=F, append=F)


########################################### DESeq2 for STAR count data #####################################

# Import gene counts table
# - skip first row (general command info)
# - make row names the gene identifiers
countdata <- read.table("./results/all_samples_raw_read_counts.txt", header = TRUE, row.names = 1)

# Remove .bam + '..' from column identifiers
colnames(countdata) <- gsub("_star.bam", "", colnames(countdata), fixed = T)

# Add gene full name for BEAVR
#countdata$ensembl <- mapIds(x = org.Mm.eg.db,keys = row.names(countdata),column = "ENSEMBL",keytype = "ENTREZID",multiVals = "first")
#remove NAs
#countdata <- na.omit(countdata)
#write to file for BEAVR
#write.table(x = as.data.frame(countdata), file = 'BEAVR_counts.txt', sep = '\t', quote = F,col.names = NA)

# Make sure ID's are correct
head(countdata)

# Add ENSEMBL for EPIC tool analysis input data formatting 
countdata$ensembl <- mapIds(x = org.Mm.eg.db,
                          keys = row.names(countdata),
                          column = "ENSEMBL",
                          keytype = "ENTREZID",
                          multiVals = "first")
countdata <- na.omit(countdata)

write.table(x = as.data.frame(countdata, normalized = T), 
            file = './EPIC/input_raw_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

#read the resulting file and write again as a csv file in order to input into jar 
EPIC_counts <- read.delim("./EPIC/input_raw_counts.txt")
write.table(x = as.data.frame(EPIC_counts), 
            file = './EPIC/csv_input_raw_counts.csv', 
            sep = ',', 
            quote = F,
            row.names = FALSE)

#read jar results and process for EÜIC
Jar_counts <- read.csv("./EPIC/jar_results.csv")

Jar_counts$GENE_NAME <- mapIds(x = org.Mm.eg.db,
                             keys = Jar_counts$ENSEMBL,
                             column = "SYMBOL",
                             keytype = "ENSEMBL",
                             multiVals = "first")

#make it all upper case
Jar_counts$GENE_NAME <- toupper(Jar_counts$GENE_NAME)

#reorder before writing to file
current_col_order <- c(colnames(Jar_counts))
col_order <- c("GENE_NAME", "WT_0h_rep_1", "WT_0h_rep_2", "WT_0h_rep_3", "WT_0h_rep_4", "WT_1h_rep_1",
               "WT_1h_rep_2", "WT_1h_rep_3", "WT_1h_rep_4", "WT_2h_rep_1", "WT_2h_rep_2", "WT_2h_rep_3",
               "WT_2h_rep_4", "WT_4h_rep_1", "WT_4h_rep_2", "WT_4h_rep_3", "WT_4h_rep_4", "ENSEMBL")
Jar_counts <- Jar_counts[, col_order]
head(Jar_counts)

#delete ensembl columns from data
Jar_counts = subset(Jar_counts, select = -c(ENSEMBL) )

#write file to .tsv
write.table(x = as.data.frame(Jar_counts), 
            file = './EPIC/Rprocessed_jar_result_raw_counts.tsv', 
            sep = '\t', 
            quote = F,
            row.names = FALSE)

# Add ENSEMBL for EPIC tool analysis input data formatting 
countdata$gene_name <- mapIds(x = org.Mm.eg.db,
                            keys = row.names(countdata),
                            column = "SYMBOL",
                            keytype = "ENTREZID",
                            multiVals = "first")
countdata <- na.omit(countdata)

write.table(x = as.data.frame(countdata, normalized = T), 
            file = './CIBERSORT/D5_mix.tsv', 
            sep = '\t', 
            quote = F,
            col.names = NA)


#Write genes > 0 (expressed) to .txt file
countdata <- countdata[apply(countdata[,-1], 1, function(x) !all(x==0)),]

write.table(x = as.data.frame(countdata, normalized = T), 
            file = 'genes_greater_zero_list.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)


# Import metadata file
# - make row names the matching sampleID's from the countdata
metadata <- read.delim("./sample.txt", row.names = 1)

# Add sampleID's to the mapping file
metadata$sampleid <- row.names(metadata)

#Add replicate column to mapping file
metadata$replicate <- c("Replicate_1", "Replicate_2", "Replicate_3", "Replicate_4", "Replicate_1", "Replicate_2", "Replicate_3", "Replicate_4","Replicate_1", "Replicate_2", "Replicate_3", "Replicate_4","Replicate_1", "Replicate_2", "Replicate_3", "Replicate_4")

# Reorder sampleID's to match featureCounts column order. 
metadata <- metadata[match(colnames(countdata), metadata$sampleid), ]

#Add replicate column to mapping file
metadata$timepoint <- c("0h", "0h", "0h", "0h", "1h", "1h", "1h", "1h","2h", "2h", "2h", "2h","4h", "4h", "4h", "4h")


# Make sure ID's are correct
head(metadata)

#Make DESeq object
# - countData : count dataframe
# - colData : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = metadata,
                                 design = ~condition)


# Find differential expressed genes
ddsMat <- DESeq(ddsMat)
resultsNames(ddsMat)


#Get basic statisics about the number of significant genes
# Get results from testing with FDR adjust pvalues for different comparisons
# generate results table for A vs ctl
res_condition_WT_1h_vs_WT_0h <- results(ddsMat, name="condition_WT_1h_vs_WT_0h", pAdjustMethod = "fdr", alpha = 0.05)
res_condition_WT_1h_vs_WT_0h

res_condition_WT_2h_vs_WT_0h <- results(ddsMat, name="condition_WT_2h_vs_WT_0h", pAdjustMethod = "fdr", alpha = 0.05)
res_condition_WT_2h_vs_WT_0h

res_condition_WT_4h_vs_WT_0h <- results(ddsMat, name="condition_WT_4h_vs_WT_0h", pAdjustMethod = "fdr", alpha = 0.05)
res_condition_WT_4h_vs_WT_0h


# Get results from testing with FDR adjust pvalues
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.01)
results <- na.omit(results)

#find number of active genes
nrow(results) #-> 14242

# Generate summary of testing. 
summary(results)


# Check directionality of the log2 fold changes
mcols(res_condition_WT_1h_vs_WT_0h, use.names = T)
mcols(res_condition_WT_2h_vs_WT_0h, use.names = T)
mcols(res_condition_WT_4h_vs_WT_0h, use.names = T)
mcols(results, use.names = T)


#List DB
edb <- org.Mm.eg.db

## List all supported keytypes.
keytypes(edb)

## List all supported columns for the select and mapIds methods.
columns(edb)

# Add gene full name
results$genesymbol <- mapIds(x = org.Mm.eg.db,
                             keys = row.names(results),
                             column = "SYMBOL",
                             keytype = "ENTREZID",
                             multiVals = "first")

results$description <- mapIds(x = org.Mm.eg.db,
                              keys = row.names(results),
                              column = "GENENAME",
                              keytype = "ENTREZID",
                              multiVals = "first")

# Add ENSEMBL
results$ensembl <- mapIds(x = org.Mm.eg.db,
                          keys = row.names(results),
                          column = "ENSEMBL",
                          keytype = "ENTREZID",
                          multiVals = "first")

# Add entrez id
results$entrez <- row.names(results)



# Subset for only significant genes (q < 0.01)
results_sig <- subset(results, padj < 0.01)
results_sig <- na.omit(results_sig)
head(results_sig)

nrow(results_sig) # -> 6742


#Write all the important results to .txt files
# Write normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(ddsMat), normalized = T), 
            file = 'normalized_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant normalized gene counts to a .txt file
write.table(x = counts(ddsMat[row.names(results_sig)], normalized = T), 
            file = 'normalized_counts_significant.txt', 
            sep = '\t', 
            quote = F, 
            col.names = NA)

# Write the annotated results table to a .txt file
write.table(x = as.data.frame(results), 
            file = "results_gene_annotated.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(results_sig), 
            file = "results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)


#Store top 3000 significant genes for DAVID analysis and write to .txt file
newdata <- results_sig[order(results_sig$log2FoldChange),]
newdata
newdata[1:3000,]

write.table(x = as.data.frame(row.names(newdata[1:3000,])), 
            file = "DAVID_results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)


################################## Stats - Exploring STAR count count data ##############################################

#explore the varying number of fragments (pairs of reads) that have been assigned to the genes for each sample
cs <- colSums(assay(ddsMat, "counts"))
hist(cs/1e6, col="grey", border="white",
     main="", xlab="column sums (per million)")

#examining the proportion of the total count for each gene
cts <- assay(ddsMat, "counts")[,c(1,4)]
idx <- rowSums(cts >= 5) == 2
cts <- cts[idx,]
p <- sweep(cts, 2, colSums(cts), "/")
mean.log.p <- rowMeans(log10(p))
hist(mean.log.p,  col="grey", border="white",
     main="", xlab="log10 proportion (geometric mean)")

#histogram
lfc <- log2(p[,2]) - log2(p[,1])
hist(lfc[between(mean.log.p, -6, -4)], 
     breaks=seq(-10,10,by=.1),
     xlim=c(-5,5),
     col="grey50", border="white",
     main="", xlab="log2 fold change of proportion")

#boxplot log2
logcounts <- log2(assay(ddsMat, "counts") + 0.1)
boxplot(logcounts, 
        main="Distribution of log counts",
        xlab="",
        ylab="Log2(raw counts+0.1)",
        las=2,cex.axis=0.8)
legend("topright", inset=c(-0.2,0), cex = 0.8,
       legend = levels(ddsMat$sampleid))


## The affy library has a density plotting function
library(affy)

## Create a list of 4 colors to use which are the same used throughout this chapter 
library(scales)
myColors <- hue_pal()(4)

par(mar = c(5, 4, 4, 8), xpd = TRUE)
## Plot the log2-transformed data with a 0.1 pseudocount
plotDensity(log2(countdata + 0.1), col=rep(myColors, each=4),
            lty=c(1:ncol(countdata)), xlab='Log2(count)',
            main='Expression Distribution')
#abline(v=-1.5, lwd=1, col='red', lty=2)

## Add a legend and vertical line
legend("topright", inset = c(- 0.82, 0),names(countdata), lty=c(1:ncol(countdata)),
       col=rep(myColors, each=4))


################################## Plotting Gene Expression Data for STAR count data ############################

#MA plot
plotMA(results, ylim = c(-5, 5))

#Plot dispersions
plotDispEsts(ddsMat)

#p-value histogram
histogram <- ggplot(as(results, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
histogram

#Plot distribution of log fold change of significant results
hist(results_sig$log2FoldChange)

#PCA
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Plot PCA by column variable only top 500 (for all the result is similar)
plotPCA(ddsMat_rlog, intgroup = "condition", ntop = 500) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 4)  + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 500 most variable genes") 


#Single gene expression plot
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Get gene with highest expression -> Cdca7 cell cycle cell division cycle associated 7"
top_gene <- rownames(results)[which.min(results$log2FoldChange)]
results["66953",]

# Plot single gene
plotCounts(dds = ddsMat, 
           gene = top_gene, 
           intgroup = "condition", 
           normalized = T, 
           transform = T)


#Heatmap
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Gather 10 significant genes and make matrix
mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:10, ]

#Find its gene symbol equivalents
collist <- row.names(mat)
#results_sig["70675", c("genesymbol")]
row.names(mat) <- c("Mrpl15", "Tcea1", "Atp6v1h" , "St18", "Pcmtd1",  "Rrs1" , "Adhfe1", "Vxn", "Mybl1", "Vcpip1")

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Condition = factor(colData(ddsMat_rlog)$condition), 
  Replicate = factor(colData(ddsMat_rlog)$replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Condition = c(WT_0h = "lightblue", WT_1h = "orange", WT_2h = "green", WT_4h = "yellow"),
  Replicate = c(Replicate_1 = "darkred", Replicate_2 = "forestgreen", Replicate_3 = "grey", Replicate_4 = "purple")
)

# Make Heatmap with pheatmap function.
## See more in documentation for customization
pheatmap(mat = mat,  
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 6.5, # Make fonts smaller
         cellwidth = 20, # Make the cells wider
         show_colnames = F, cluster_cols=FALSE)

#Plot top 10 genes from BEAVR-result in heatmap
gene_symbol_annot <- c("Cmpk2", "Mx1", "Rsad2", "Nos2", "Saa3", "Ccl5", "Acod1", "Cxcl10", "Il1b", "Ptgs2") #gene symbols
marker.genes <- c("22169", "17857", "58185", "18126", "20210", "20304", "16365", "15945", "16176", "19225") #entrez ids

# Gather 10 significant genes and make matrix
mat <- assay(ddsMat_rlog[marker.genes])

#Set rownames of mat to gene symbols
row.names(mat) <- gene_symbol_annot

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Timepoint = factor(colData(ddsMat_rlog)$timepoint), 
  Replicate = factor(colData(ddsMat_rlog)$replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Timepoint = c("0h" = "lightblue", "1h" = "orange", "2h" = "green", "4h" = "yellow"),
  Replicate = c(Replicate_1 = "darkred", Replicate_2 = "forestgreen", Replicate_3 = "grey", Replicate_4 = "purple")
)

#Heatmap for various (atherosclerosis) biomarkers showing time course progression
pheatmap(mat = mat,  
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 7.5, # Make fonts smaller
         cellwidth = 15, # Make the cells wider
         show_colnames = F, cluster_cols=FALSE)


#Volcanoplot
# Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
## - Change pvalues to -log10 (1.3 = 0.05)
data <- data.frame(gene = row.names(results),
                   pval = -log10(results$padj), 
                   lfc = results$log2FoldChange)

# Remove any rows that have NA as an entry
data <- na.omit(data)

# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < 0 & data$pval > 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "nonsignificant"))

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

# Add ggplot2 layers
vol +   
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "red", Decreased = "blue", nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("conditions"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values


################################## Finding Pathways from Differential Expressed Genes for STAR counts ############################

# Remove any genes that do not have any entrez identifiers
results_sig_entrez <- subset(results_sig, is.na(entrez) == FALSE)
na.omit(results_sig_entrez)

# Create a matrix of gene log2 fold changes
gene_matrix <- results_sig_entrez$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(gene_matrix) <- results_sig_entrez$entrez

# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(gene_matrix)


#Enrich genes using the KEGG database
kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'mouse',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)
head(kegg_enrich)

# Plot results
barplot(kegg_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)


#Gene ontology
go_enrich <- enrichGO(gene = names(gene_matrix),
                      OrgDb = 'org.Mm.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.01, 
                      qvalueCutoff = 0.10)

head(go_enrich)

#write GO results to file
write.table(x = as.data.frame(go_enrich), 
            file = "./Resulting_tables/GO_enrich_result.txt", 
            sep = '\t', 
            quote = F)

#dotplot
dotplot(go_enrich, showCategory = 10, font.size=14)


# Plot results
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)


#Plotting KEGG pathways
# Plot specific KEGG pathways (with fold change) 
## pathway.id : KEGG pathway identifier
pathview(gene.data = gene_matrix, 
         pathway.id = "04070", 
         species = "mouse")



#Active pathways for Signaling Pathway Impact Analysis (>20 DEGs)

# Get a vector of log(FC) values for all significant genes
sig_genes <- subset(results, padj < 0.05, select = log2FoldChange)[[1]]
# Make it a named vector by assigning the Entrez ID's to each log(FC) value
names(sig_genes) <- subset(results, padj < 0.05, select = entrez)[[1]]
# A complete list of Entrez IDs for all genes in this experiment
all_genes <- results$entrez

pander(all_genes[which(idx)[100:110]])

# Process all signaling pathways to see if they are inhibited or activated
spia_result <- spia(de=sig_genes, all=all_genes, organism="mmu", plots=TRUE)
head(spia_result)

#write to .txt
write.table(x = as.data.frame(spia_result), 
            file = "./SPIA/spia_result.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

#plot
plotP(spia_result)


#Read spia_result table
spia_result <- read.delim("./SPIA/spia_result.txt", header = TRUE, row.names = 1)



#Prepare data for import to iPathwayGuide tool
prep_ipathway <- results_sig
prep_ipathway <- na.omit(prep_ipathway)

#make rownames from entrez id to gene symbol
prep_ipathway$genesymbol <- mapIds(x = org.Mm.eg.db,
                             keys = row.names(prep_ipathway),
                             column = "SYMBOL",
                             keytype = "ENTREZID",
                             multiVals = "first")

#remove missing values 
prep_ipathway <- na.omit(prep_ipathway)

#change rownames to gene symbol
rownames(prep_ipathway) <- prep_ipathway$genesymbol

#remove genesymbol column
prep_ipathway = subset(prep_ipathway, select = -c(genesymbol))

#write output to .csv file
write.table(x = as.data.frame(prep_ipathway), 
            file = "./iPathwayGuide/Deseq_res.csv", 
            sep = ',', 
            quote = F, col.names = NA)


######################################## Mapping reads to chromosomes for STAR counts ###################################

#Did the same for all bam.bai files(!)
output.bam.index <- "./all_samples_STAR_bam_bai/WT_0h_rep_3_star.bam.bai"
output.bamFile <- "./data_all_samples_STAR/WT_0h_rep_3_star.bam"


#Find number of mapped reads per chromosome
chr.mapping.stats <- idxstatsBam(output.bamFile, index=output.bam.index)
chr.mapping.stats

#Barplot percent reads
total.mapped.reads <- sum(chr.mapping.stats$mapped)
chr.mapping.stats$mapped.prop <- chr.mapping.stats$mapped/total.mapped.reads* 10
barplot(chr.mapping.stats$mapped.prop, 
        names.arg=as.character(chr.mapping.stats$seqnames), las=2, col='slateblue', ylim=c(0,2)) + title("Mapped reads per chromosome")

#Barplot
rownames(chr.mapping.stats) <- chr.mapping.stats$seqnames
barplot(chr.mapping.stats$mapped, 
        names.arg=as.character(chr.mapping.stats$seqnames), col='slateblue')

################################# Mapper comparison plots ###################################

#### hisat ####

# Import gene counts table
countdata_hisat <- read.table("./results/hisat/hisat_all_samples_raw_read_counts.txt", header = TRUE, row.names = 1)

# Remove .bam + '..' from column identifiers
colnames(countdata_hisat) <- gsub("_hisat.bam", "", colnames(countdata_hisat), fixed = T)

# Make sure ID's are correct
head(countdata_hisat)

#Make DESeq object
ddsMat_hisat <- DESeqDataSetFromMatrix(countData = countdata_hisat,
                                       colData = metadata,
                                       design = ~condition)


# Find differential expressed genes
ddsMat_hisat <- DESeq(ddsMat_hisat)
resultsNames(ddsMat_hisat)

# Get results from testing with FDR adjust pvalues
results_hisat <- results(ddsMat_hisat, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
summary(results_hisat)

#MA plot
plotMA(results_hisat, ylim = c(-5, 5))

#Plot dispersions
plotDispEsts(ddsMat_hisat)

#p-value histogram
histogram <- ggplot(as(results_hisat, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.05, fill = "Royalblue", boundary = 0)
histogram



#### contextmap ####

# Import gene counts table
countdata_contextmap <- read.table("./results/contextmap/contextmap_all_samples_raw_read_counts.txt", header = TRUE, row.names = 1)

# Remove .bam + '..' from column identifiers
colnames(countdata_contextmap) <- gsub("_contextmap.bam", "", colnames(countdata_contextmap), fixed = T)

# Make sure ID's are correct
head(countdata_contextmap)

#Make DESeq object
ddsMat_contextmap <- DESeqDataSetFromMatrix(countData = countdata_contextmap,
                                            colData = metadata,
                                            design = ~condition)


# Find differential expressed genes
ddsMat_contextmap <- DESeq(ddsMat_contextmap)
resultsNames(ddsMat_contextmap)

# Get results from testing with FDR adjust pvalues
results_contextmap <- results(ddsMat_contextmap, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
summary(results_contextmap)

#MA plot
plotMA(results_contextmap, ylim = c(-5, 5))

#Plot dispersions
plotDispEsts(ddsMat_contextmap)

#p-value histogram
histogram <- ggplot(as(results_contextmap, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.05, fill = "Royalblue", boundary = 0)
histogram


#### tophat2 ####

# Import gene counts table
countdata_tophat2 <- read.table("./results/tophat2/tophat2_all_samples_raw_read_counts.txt", header = TRUE, row.names = 1)

# Remove .bam + '..' from column identifiers
colnames(countdata_tophat2) <- gsub("_tophat2.bam", "", colnames(countdata_tophat2), fixed = T)

# Make sure ID's are correct
head(countdata_tophat2)

#Make DESeq object
ddsMat_tophat2 <- DESeqDataSetFromMatrix(countData = countdata_tophat2,
                                         colData = metadata,
                                         design = ~condition)


# Find differential expressed genes
ddsMat_tophat2 <- DESeq(ddsMat_tophat2)
resultsNames(ddsMat_tophat2)

# Get results from testing with FDR adjust pvalues
results_tophat2 <- results(ddsMat_tophat2, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
summary(results_tophat2)

#MA plot
plotMA(results_tophat2, ylim = c(-5, 5))

#Plot dispersions
plotDispEsts(ddsMat_tophat2)

#p-value histogram
histogram <- ggplot(as(results_tophat2, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.05, fill = "Royalblue", boundary = 0)
histogram

################################### EdgeR and DESeq comparison (Venn diagram) ################################

rawCountTable <- read.table("./results/all_samples_raw_read_counts.txt", header=TRUE, sep="\t", row.names=1)
colnames(rawCountTable) <- gsub("_star.bam", "", colnames(rawCountTable), fixed = T)
sampleInfo <- read.table("./sample.txt", header=TRUE, sep="\t", row.names=1)

dgeFull <- DGEList(rawCountTable, group=sampleInfo$condition)
dgeFull

#remove genes with zero counts for all samples
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
group=dgeFull$samples$group)
head(dgeFull$counts)

#normalizing
dgeFull <- calcNormFactors(dgeFull, method="TMM")
#dgeFull$samples

dgeFull <- estimateDisp(dgeFull)
sqrt(dgeFull$common.dispersion) # biological coefficient of variation
#plotBCV(dgeFull)

dgeFull$samples$group <- c("WT_0h","WT_0h", "WT_0h", "WT_0h", "WT_1h", "WT_1h", "WT_1h", "WT_1h", "WT_2h", "WT_2h", "WT_2h",
                       "WT_2h", "WT_4h", "WT_4h", "WT_4h", "WT_4h")
et <- exactTest(dgeFull)
results_edgeR <- topTags(et, n = nrow(rawCountTable), sort.by = "none")
head(results_edgeR$table)

#volcanoplot
volcanoData <- cbind(results_edgeR$table$logFC, -log10(results_edgeR$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)
plot(volcanoData, pch=19)

# Filter on adjusted p-value and get the rownames
# edgeR
edgeR.degs <- rownames(results_edgeR)[results_edgeR$table$FDR < 0.01]
sum(results_edgeR$table$FDR < 0.01) # -> 1595

#plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < 0.05])
#abline(h = c(-2, 2), col = "blue")

# DESeq2
DESeq.degs <- row.names(results_sig)
sum(results_sig$padj < 0.01)



# Venn-diagram 
# Calculate the intersection of the two sets
deg.intersect = length(intersect(edgeR.degs, DESeq.degs))
deg.venn <- list('intersect' = deg.intersect,
                 'edgeR' = length(edgeR.degs),
                 'DESeq2' = length(DESeq.degs))

# Arguments for a pairwise (two-sets) venn-diagram are sizes for set1, set2 and overlap (intersect)
# Many more functions are available for triple, quad and quantuple diagrams (starting with 'draw.***')
venn.plot <- draw.pairwise.venn(deg.venn$edgeR, deg.venn$DESeq2, deg.venn$intersect,
                                category = c("edgeR", "DESeq2"), scaled = F,
                                fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                                cat.pos = c(0, 0))

# Actually plot the plot
grid.draw(venn.plot)

#boxplot of log-fold changes
d1 <- results_edgeR$table$logFC
d2 <- results_sig$log2FoldChange
boxplot(d1, main="log-FC values for edgeR",
xlab="edgeR", ylab="log FC")
boxplot(d2, main="log-FC values for DESeq2",
        xlab="DESeq2", ylab="log FC")


#Combined boxplot
boxplot(results_edgeR$table$logFC[results_edgeR$table$FDR < 0.01], 
        results_sig$log2FoldChange[results_sig$padj <= 0.01],
        names=c('edgeR', 'DESeq2'), ylab="log-FC",
        main = 'log-FC values for all DEGs by package')
abline(h=0, lty=2, col='red')
###################################### Biomarker analysis #####################################

#Read results_sig with gene symbol as row name
restructured_results_sig <- read.delim("./restructured_results_sig.txt", row.names = 1)

# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

###List of marker genes: (total of 19 genes)
##Chemokines: Ccl3, Ccl4 (CC subfamily) Cxcl10 (CXC subfamily)
##Pro-inflammatory cytokines: Il1b  Il6, Il6st (Hematopoesis) Tnfaip8l2, Tnfaip6, (TNF family)
##Anti-inflammatory cytokines: Il10, Il10ra (IL-10 family) Tgfbr3, Tgfbr1, Tgfbrap1 (TGF beta family)
##Atherosclerosis specifics: Ifng (Interferon gamma), Il2ra (interleukin 2), Csf2 (Granulocyte Macrophage Colony Stimulating Factor (GMCSF)), Il18, Il12a, Il12b

gene_symbol_annot <- c("Ccl3", "Ccl4", "Cxcl10", "Il1b", "Il6", "Il6st", "Tnfaip8l2", "Tnfaip6", "Il10", "Il10ra", "Tgfbr3", "Tgfbr1", "Tgfbrap1", "Ifng", "Csf2", "Il2ra", "Il18", "Il12a", "Il12b") #gene symbols
marker.genes <- c("20302", "20303", "15945", "16176", "16193", "16195", "69769", "21930", "16153", "16154", "21814", "21812", "73122", "15978", "12981", "16184", "16173", "16159", "16160") #entrez ids

# Gather 10 significant genes and make matrix
mat <- assay(ddsMat_rlog[marker.genes])

#Set rownames of mat to gene symbols
row.names(mat) <- gene_symbol_annot

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Timepoint = factor(colData(ddsMat_rlog)$timepoint), 
  Replicate = factor(colData(ddsMat_rlog)$replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Timepoint = c("0h" = "lightblue", "1h" = "orange", "2h" = "green", "4h" = "yellow"),
  Replicate = c(Replicate_1 = "darkred", Replicate_2 = "forestgreen", Replicate_3 = "grey", Replicate_4 = "purple")
)

#Heatmap for various (atherosclerosis) biomarkers showing time course progression
pheatmap(mat = mat,  
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 7.5, # Make fonts smaller
         cellwidth = 15, # Make the cells wider
         show_colnames = F, cluster_cols=FALSE, cluster_rows = FALSE)

##################################### Biomarker part 2 #####################
###List of marker genes: 
gene_symbol_annot <- c("Ccl3", "Ccl4", "Ccl5", "Cxcl10", "Il1b", "Vegfa", "Tnfa", "Ifng", "Il12", "Il18", "Il2", "Il10", "Tgfb", "Il4", "Il13", "Il6") #gene symbols
marker.genes <- c("20302", "20303", "20304", "15945", "16176", "22339", "21928", "15978", "16160", "16173", "16186", "16153", "73122", "16190", "16164", "16193") #entrez ids

# Gather 10 significant genes and make matrix
mat <- assay(ddsMat_rlog[marker.genes])

#Set rownames of mat to gene symbols
row.names(mat) <- gene_symbol_annot

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Timepoint = factor(colData(ddsMat_rlog)$timepoint), 
  Replicate = factor(colData(ddsMat_rlog)$replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Timepoint = c("0h" = "lightblue", "1h" = "orange", "2h" = "green", "4h" = "yellow"),
  Replicate = c(Replicate_1 = "darkred", Replicate_2 = "forestgreen", Replicate_3 = "grey", Replicate_4 = "purple")
)

#Heatmap for various (atherosclerosis) biomarkers showing time course progression
pheatmap(mat = mat,  
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 7.5, # Make fonts smaller
         cellwidth = 15, # Make the cells wider
         show_colnames = F, cluster_cols=FALSE, cluster_rows = FALSE)

####################################### Time course analysis #####################################

#find out what is the gene with smallest p-adj value
smallest_padj_gene <- rownames(results)[which.min(results$padj)] 
results["12633", ] #Cflar: "CASP8 and FADD-like apoptosis regulator" 

#check if this is a unique value or not
results[which(results$padj == min(results$padj)) ] # -> this is not unique

#plot counts
time.p.counts <- plotCounts(ddsMat, which.min(results$padj), 
                   intgroup = c("condition","timepoint"), returnData = TRUE)

ggplot(time.p.counts,
       aes(x = timepoint, y = count, color = condition, group = condition)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10() + ggtitle("Normalized counts for the gene with smallest p-adj value (Cflar) to show condition-specific changes over time")


####################################### Biotype of detected genes #################################

#read the biotypes file in and map it to the counts data column

#Read results_sig with gene symbol as row name
biotypes_id_list <- read.table("./EPIC/mahima_geneBiotypes.tsv", header = TRUE)

#add entrz id
biotypes_id_list$entrez <- mapIds(x = org.Mm.eg.db,
                             keys = biotypes_id_list$gene_id,
                             column = "ENTREZID",
                             keytype = "ENSEMBL",
                             multiVals = "first")

write.table(x = as.data.frame(biotypes_id_list), 
            file = "blah.list", 
            sep = '\t', 
            quote = F,
            col.names = NA)

#read the resulting biotypes
biotypes_id_list <- read.table("./biotypes_with_entrez_id.txt", header = TRUE)


#map the biotypes to ddsMat count data to plot it
ddsMat_biotype <- assay(ddsMat, "counts") 
ddsMat_biotype <- as.data.frame(ddsMat_biotype, normalized = T)
biotyple_count_table <- merge(ddsMat_biotype, biotypes_id_list, by = 0, all = FALSE)

#melt the table
library(data.table)
biotyple_count_table <- reshape2::melt(biotyple_count_table, id.vars = c("Row.names", "biotype"), 
                             measure.vars = c("WT_0h_rep_1", "WT_0h_rep_2", 
                                              "WT_0h_rep_3", "WT_0h_rep_4", "WT_1h_rep_1",
                                              "WT_1h_rep_2", "WT_1h_rep_3", "WT_1h_rep_4", 
                                              "WT_2h_rep_1", "WT_2h_rep_2",
                                              "WT_2h_rep_3", "WT_2h_rep_4", "WT_4h_rep_1", 
                                              "WT_4h_rep_2", "WT_4h_rep_3", "WT_4h_rep_4"),
     variable.name = "sample", value.name = "counts")

#remove unneccessary biotypes
biotyple_count_table_final <- biotyple_count_table[!(biotyple_count_table$biotype == "IG_C_gene" | 
                                                    biotyple_count_table$biotype =="IG_D_gene" |
                                                    biotyple_count_table$biotype =="IG_J_gene" |
                                                    biotyple_count_table$biotype =="IG_V_gene" |
                                                    biotyple_count_table$biotype =="IG_V_pseudogene" |
                                                    biotyple_count_table$biotype =="misc_RNA" |
                                                    biotyple_count_table$biotype =="polymorphic_pseudogene" |
                                                    biotyple_count_table$biotype =="pseudogene" |
                                                    biotyple_count_table$biotype =="ribozyme" |
                                                    biotyple_count_table$biotype =="scaRNA" |
                                                    biotyple_count_table$biotype =="snoRNA" |
                                                    biotyple_count_table$biotype =="TR_C_gene" |
                                                    biotyple_count_table$biotype =="TR_D_gene" |
                                                    biotyple_count_table$biotype =="TR_J_gene" |
                                                    biotyple_count_table$biotype =="TR_V_gene" |
                                                    biotyple_count_table$biotype =="TR_V_pseudogene" |
                                                    biotyple_count_table$biotype =="transcribed_unitary_pseudogene" |
                                                    biotyple_count_table$biotype =="transcribed_unprocessed_pseudogene" |
                                                    biotyple_count_table$biotype =="translated_unprocessed_pseudogene"),]


#Plot stacked barplot
bt <- ggplot(biotyple_count_table_final, aes(fill=biotype, y=counts, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) 

bt

#barplot of percentage
ggplot(biotyple_count_table_final, aes(fill=biotype, y= (counts/sum(counts)) * 10, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) + xlab("samples") +
  ylab("count in percent")



####################################### Generate all table outputs ################################

#Read  D5 mix
D5_count <- read.table("./CIBERSORT/D5_mix.tsv", header = TRUE, row.names = 1)


#Write genes > 0 (expressed) to .txt file
D5_count <- D5_count[apply(D5_count[,-1], 1, function(x) !all(x==0)),]

#nrow
nrow(D5_count) #20336 genes


#Make deseq object
ddsMat <- DESeqDataSetFromMatrix(countData = D5_count,
                                 colData = metadata,
                                 design = ~condition)


# Find differential expressed genes
ddsMat <- DESeq(ddsMat)


# Get results from testing with FDR adjust pvalues
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)
# Add ENSEMBL
results$entrez <- mapIds(x = org.Mm.eg.db,
                          keys = row.names(results),
                          column = "ENTREZID",
                          keytype = "SYMBOL",
                          multiVals = "first")
results <- na.omit(results)


# Generate summary of testing. 
summary(results)

#now 
nrow(results) #15604


#View
View(as.data.frame(results))


#Write all expressed genes
write.table(x = as.data.frame(results), 
            file = "./Resulting_tables/Set_of_expressed_genes.tsv", 
            sep = '\t', 
            quote = F,
            col.names = NA)



# Subset for only significant genes (q < 0.05)
results_sig <- subset(results, padj < 0.05)
results_sig <- na.omit(results_sig)
head(results_sig)
nrow(results_sig) #8168

#Write all differentially expressed genes
write.table(x = as.data.frame(results_sig), 
            file = "./Resulting_tables/Set_of_differentially_expressed_genes.tsv", 
            sep = '\t', 
            quote = F,
            col.names = NA)

#View
View(as.data.frame(results_sig))


#Read alternative splicing results
dexSeq_results <- read.delim("./Resulting_tables/dexSeq_results.tsv", header = TRUE)

#View
View(dexSeq_results) #360917

#nrow
nrow(dexSeq_results)


#Read differential gene sets
GO_enrich_result <- read.delim("./Resulting_tables/GO_enrich_result.txt", header = TRUE, row.names = 1)

#View
View(GO_enrich_result)

#nrow
nrow(GO_enrich_result) #3794


#Aktive genesets
# Create a matrix of gene log2 fold changes
gene_matrix <- results$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(gene_matrix) <- results$entrez

# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(gene_matrix)


#Enrich genes using the KEGG database
kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'mouse',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)
head(kegg_enrich)

# Plot results
barplot(kegg_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)


#Gene ontology
go_enrich <- enrichGO(gene = names(gene_matrix),
                      OrgDb = 'org.Mm.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

head(go_enrich)

#write GO results to file
write.table(x = as.data.frame(go_enrich), 
            file = "./Resulting_tables/GO_enrich_active_genesets_result.txt", 
            sep = '\t', 
            quote = F)

GO_enrich_active_genes_result <- read.delim("./Resulting_tables/GO_enrich_active_genesets_result.txt", header = TRUE, row.names = 1)
View(GO_enrich_active_genes_result)
nrow(GO_enrich_active_genes_result) #4638

#dotplot
dotplot(go_enrich, showCategory = 10, font.size=14)


# Plot results
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)


#Active pathways for Signaling Pathway Impact Analysis (>20 DEGs)
# A complete list of Entrez IDs for all genes in this experiment
all_genes <- results$entrez

# Process all signaling pathways to see if they are inhibited or activated
spia_result <- spia(all=all_genes, organism="mmu", plots=TRUE)
head(spia_result)

#write to .txt
write.table(x = as.data.frame(spia_result), 
            file = "./SPIA/spia_result.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

#plot
plotP(spia_result)


#Read spia_result table
spia_result <- read.delim("./SPIA/spia_result.txt", header = TRUE, row.names = 1)


################################ R-session info ###########################################

#Write R session info to file 
writeLines(capture.output(sessionInfo()), "C:/Users/mahim/GoBi/Blockteil/Paper/Supplementary_files/sessionInfo.txt")

