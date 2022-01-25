#read file
Corona_CountTable = read.table("../gene.counts.star", header=TRUE, row.names=1 )
head( Corona_CountTable )
counts = as.matrix(Corona_CountTable) #as matrix
dim(counts)

#remove .bam/.sam extension from headers
colnames(Corona_CountTable) <- gsub("\\.[sb]am$", "", colnames(Corona_CountTable))
#colnames(Corona_CountTable)

#meta data
Corona_MetaData <- read.table("sample.list", header=TRUE)
head(Corona_MetaData)
#colnames(Corona_MetaData)
#condition = factor( Corona_MetaData$condition)
#condition


#DeSeq
library( "DESeq" )
cds = newCountDataSet( countTable, condition )












#############################################################################

#############################################################################
#load libraries
#BiocManager::install(c("DESeq2", "clusterProfiler", "biomaRt", "ReactomePA","DOSE", "KEGG.db", "org.Hs.eg.db", "genefilter", "GO.db","topGO", "gage", "ggsci"))

library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
library(KEGG.db)
library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)


# Import gene counts table
# - skip first row (general command info)
# - make row names the gene identifiers
Corona_CountTable = read.table("../gene.counts.star", header=TRUE, row.names=1 )
colnames(Corona_CountTable) <- gsub("\\.[sb]am$", "", colnames(Corona_CountTable))
head( Corona_CountTable )

# Import metadata file
# - make row names the matching sampleID's from the countdata
metadata <- read.delim("sample.list", row.names = 1)

# Add sampleID's to the mapping file
metadata$sampleid <- row.names(metadata)

# Make sure ID's are correct
head(metadata)


# - countData : count dataframe
# - colData : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = metadata,
                                 design = ~condition)


# Find differential expressed genes
ddsMat <- DESeq(ddsMat)


# Get results from testing with FDR adjust pvalues
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
summary(results)


# Check directionality of the log2 fold changes
## Log2 fold change is set as (LoGlu / HiGlu)
## Postive fold changes = Increased in LoGlu
## Negative fold changes = Decreased in LoGlu
mcols(results, use.names = T)


################################ PLOTTING GENE EXPRESSION DATA #########################
