# downloading gtex data with recount3
library(recount3)
library(DESeq2)
library(edgeR)
library(RUVSeq)
library("snapcount")
library(pheatmap)
library(org.Mm.eg.db)
library(data.table)

## get table of available projects from recount3
human_projects <- available_projects()

dim(human_projects)
head(human_projects)

## Find all the mouse projects
mouse_projects <- available_projects(organism = "mouse")

## Explore the results
dim(mouse_projects)
head(mouse_projects)

####################################################################### SRP134109 - Il1b project  #########################################################

## Find the project 
proj_info <- subset(
  mouse_projects,
  project == "SRP134109" & project_type == "data_sources"
)

## Create a RangedSummarizedExperiment (RSE) object at the gene level
rse_gene_SRP134109 <- create_rse(proj_info)

## Explore that RSE object
rse_gene_SRP134109

#save RSE object as rds file
saveRDS(rse_gene_SRP134109, file = "data_SRP134109/rse_gene_SRP134109.Rds")

rse_gene_SRP134109 <- readRDS("data_SRP134109/rse_gene_SRP134109.Rds")

#Same thing 
data_SRP134109 <- recount3::create_rse_manual(
  project = "SRP134109",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)


#inspect rse
data_SRP134109

#save RSE object as rds file
saveRDS(data_SRP134109, file = "data_SRP134109/Gene_data_SRP134109.Rds")

data_SRP134109 <- readRDS("data_SRP134109/Gene_data_SRP134109.Rds")

                                          ##################### Gene Bp-level count data ###################################

# Once you have your RSE object, you can transform the raw coverage
## base-pair coverage counts using transform_counts().
assay(rse_gene_SRP134109, "counts") <- transform_counts(rse_gene_SRP134109)


# load Metadata information 
#Sample_title	"SL139389: IgG_Control_1"	"SL139390: IgG_Control_2"	"SL139391: IL1_AB_1"	"SL139392: IL1_AB_2"	"SL139393: IgG_Control_3"	"SL139394: IgG_Control_4"	"SL139395: IL1_AB_3"	"SL139396: IL1_AB_4"
metadata <- read.delim("./data_SRP134109/metadata.txt", row.names = 1)


#Deseq2 we compute log fold changes for the condition IgG vs. Il1b treatment
ddsMat_SRP134109 <- DESeqDataSetFromMatrix(countData = assay(rse_gene_SRP134109, "counts"),
                                 colData = metadata,
                                 design = ~condition)


# Find differential expressed genes
ddsMat_SRP134109 <- DESeq(ddsMat_SRP134109)

# Get results from testing with FDR adjust pvalues
results_SRP134109 <- results(ddsMat_SRP134109, pAdjustMethod = "fdr", alpha = 0.01)

#remove NA's
results_SRP134109 <- na.omit(results_SRP134109)

# Generate summary of testing. 
summary(results_SRP134109)

#Write it to file
write.table(x = as.data.frame(results_SRP134109), 
            file = './results_SRP134109/Gene_logfoldchange_values_deseq.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)


                                          ##################### Exon Bp-level count data ###################################

data_SRP134109_exon <- recount3::create_rse_manual(
  project = "SRP134109",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "exon"
)

#inspect rse
data_SRP134109_exon
rowRanges(data_SRP134109_exon)

#save RSE object as rds file
saveRDS(data_SRP134109_exon, file = "data_SRP134109/Exon_data_SRP134109.Rds")

#TRead RDS
data_SRP134109_exon <- readRDS("data_SRP134109/Exon_data_SRP134109.Rds")

#Try scaling data
data_SRP134109_exon <- transform_counts(data_SRP134109_exon)


# load Metadata information 
#Sample_title	"SL139389: IgG_Control_1"	"SL139390: IgG_Control_2"	"SL139391: IL1_AB_1"	"SL139392: IL1_AB_2"	"SL139393: IgG_Control_3"	"SL139394: IgG_Control_4"	"SL139395: IL1_AB_3"	"SL139396: IL1_AB_4"
metadata <- read.delim("./data_SRP134109/metadata.txt", row.names = 1)



#Deseq2 we compute log fold changes for the condition IgG vs. Il1b treatment
ddsMat_SRP134109_exon <- DESeqDataSetFromMatrix(countData = data_SRP134109_exon,
                                           colData = metadata,
                                           design = ~condition)


# Find differential expressed genes
ddsMat_SRP134109_exon <- DESeq(ddsMat_SRP134109_exon)

# Get results from testing with FDR adjust pvalues
results_SRP134109_exon <- results(ddsMat_SRP134109_exon, pAdjustMethod = "fdr", alpha = 0.05)

#remove NA's
results_SRP134109_exon <- na.omit(results_SRP134109_exon)

# Generate summary of testing. 
summary(results_SRP134109_exon)

#Write it to file
write.table(x = as.data.frame(results_SRP134109_exon), 
            file = './results_SRP134109/Exon_logfoldchange_values_deseq.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)


                     ##################### Junction-level count data ###################################

data_SRP134109_jxn <- recount3::create_rse_manual(
  project = "SRP134109",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "jxn"
)

#inspect rse
data_SRP134109_jxn
rowRanges(data_SRP134109_jxn)

#save RSE object as rds file
saveRDS(data_SRP134109_jxn, file = "data_SRP134109/Junction_data_SRP134109.Rds")

#TRead RDS
data_SRP134109_jxn <- readRDS("data_SRP134109/Junction_data_SRP134109.Rds")

#Try scaling data
data_SRP134109_jxn <- transform_counts(data_SRP134109_jxn)


#Deseq2 we compute log fold changes for the condition IgG vs. Il1b treatment
ddsMat_SRP134109_jxn <- DESeqDataSetFromMatrix(countData = data_SRP134109_jxn,
                                                colData = metadata,
                                                design = ~condition)


# Find differential expressed genes
ddsMat_SRP134109_jxn <- DESeq(ddsMat_SRP134109_jxn)

# Get results from testing with FDR adjust pvalues
results_SRP134109_jxn <- results(ddsMat_SRP134109_jxn, pAdjustMethod = "fdr", alpha = 0.05)

#remove NA's
results_SRP134109_jxn <- na.omit(results_SRP134109_jxn)

#Write it to file
write.table(x = as.data.frame(results_SRP134109_jxn), 
            file = './results_SRP134109/Junction_logfoldchange_values_deseq.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)


           ##################### Basepair-level count data (BigWig sample level files) ###################################

basename(head(rse_gene_SRP134109$BigWigURL))
#rtracklayer::import.bw()

############################################# Snapcount for SRP134109 for Il1b #################################

#query-builder class
#query gene
sb <- QueryBuilder(compilation="gtex", regions="IL1B")
Il1b.gene <- query_gene(sb)
dim(Il1b.gene)
head(Il1b.gene)
assay(Il1b.gene)

#query exon
Il1b.exon <- query_exon(sb)
dim(Il1b.exon)
head(Il1b.exon)
assay(Il1b.exon)

#query junctions
Il1b.jx.all <- query_jx(sb)
dim(Il1b.jx.all)
head(Il1b.jx.all)
assay(Il1b.jx.all)


############################################ Get top X genes to investigate for DS5 #############################################

Dif_expr_genes_DS5 <- read.table("C:/Users/mahim/GoBi/Blockteil/Analysis_proj5/Resulting_tables/Set_of_differentially_expressed_genes.tsv", header = TRUE, row.names = 1)

#order padj from smallest to largest and print head
Dif_expr_genes_DS5 <-Dif_expr_genes_DS5[order(Dif_expr_genes_DS5$padj, Dif_expr_genes_DS5$pvalue),]
head(Dif_expr_genes_DS5)

#subset to see all genes with p-value and padj being 0
Dif_expr_genes_DS5[(Dif_expr_genes_DS5[,"pvalue"]==0|Dif_expr_genes_DS5[,"padj"]==0),]

#print hand-picked 5 genes from this list to compare with other recount3 datasets 
names_of_genes_of_interest <- c("Klf2", "Klf4", "Nos3")
genes_of_interest_listing <- subset(Dif_expr_genes_DS5, rownames(Dif_expr_genes_DS5) %in% names_of_genes_of_interest)
View(genes_of_interest_listing)

#make heatmap of these 5 genes for DS5
# Import gene counts table
# - skip first row (general command info)
# - make row names the gene identifiers
countdata <- read.table("C:/Users/mahim/GoBi/Blockteil/Analysis_proj5/CIBERSORT/D5_mix.tsv", header = TRUE, row.names = 1)

# Remove .bam + '..' from column identifiers
colnames(countdata) <- gsub("_star.bam", "", colnames(countdata), fixed = T)

#Write genes > 0 (expressed) to .txt file
countdata <- countdata[apply(countdata[,-1], 1, function(x) !all(x==0)),]

# Import metadata file
# - make row names the matching sampleID's from the countdata
metadata <- read.delim("C:/Users/mahim/GoBi/Blockteil/Analysis_proj5/sample.txt", row.names = 1)

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

# Get results from testing with FDR adjust pvalues
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.01)

# Subset for only significant genes (q < 0.05)
results_sig <- subset(results, padj < 0.01)
results_sig <- na.omit(results_sig)
head(results_sig)

#Heatmap
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Gather 10 significant genes and make matrix
mat <- assay(ddsMat_rlog[names_of_genes_of_interest])

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



###################### Region data - Check for these genes in SRP134109 gene-level count data #############

rse_gene_SRP134109 <- readRDS("data_SRP134109/rse_gene_SRP134109.Rds")
assay(rse_gene_SRP134109, "counts") <- transform_counts(rse_gene_SRP134109)

#Clip NKS from count data
rownames(rse_gene_SRP134109) <- gsub("\\..*","",rownames(rse_gene_SRP134109))

# load Metadata information 
#Sample_title	"SL139389: IgG_Control_1"	"SL139390: IgG_Control_2"	"SL139391: IL1_AB_1"	"SL139392: IL1_AB_2"	"SL139393: IgG_Control_3"	"SL139394: IgG_Control_4"	"SL139395: IL1_AB_3"	"SL139396: IL1_AB_4"
metadata <- read.delim("./data_SRP134109/metadata.txt", row.names = 1)


#Deseq2 we compute log fold changes for the condition IgG vs. Il1b treatment
ddsMat_SRP134109 <- DESeqDataSetFromMatrix(countData = assay(rse_gene_SRP134109, "counts"),
                                           colData = metadata,
                                           design = ~condition)


# Find differential expressed genes
ddsMat_SRP134109 <- DESeq(ddsMat_SRP134109)

# Get results from testing with FDR adjust pvalues
results_SRP134109 <- results(ddsMat_SRP134109, pAdjustMethod = "fdr", alpha = 0.01)

#remove NA's
results_SRP134109 <- na.omit(results_SRP134109)

# Subset for only significant genes (q < 0.05)
results_SRP134109_sig <- subset(results_SRP134109, padj < 0.01)
results_SRP134109_sig <- na.omit(results_SRP134109_sig)
head(results_SRP134109_sig)

#Make heatmap showing these 5 genes of interest

#Heatmap
# Convert all samples to rlog
ddsMat_rlog_SRP134109 <- rlog(ddsMat_SRP134109, blind = FALSE)

# Gather 10 significant genes and make matrix
ensembllist <- mapIds(x = org.Mm.eg.db,
                      keys = names_of_genes_of_interest,
                      column = "ENSEMBL",
                      keytype = "SYMBOL",
                      multiVals = "first")

#names_of_genes_of_interest_ensembl <- c("ENSMUSG00000027398", "ENSMUSG00000027776", "ENSMUSG00000055170", "ENSMUSG00000021281", "ENSMUSG00000004296")
names_of_genes_of_interest_ensembl <- c("ENSMUSG00000055148", "ENSMUSG00000003032", "ENSMUSG00000028978")


mat <- assay(ddsMat_rlog_SRP134109[names_of_genes_of_interest_ensembl])
row.names(mat) <- c("Klf2", "Klf4", "Nos3")

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Condition = factor(colData(ddsMat_SRP134109)$condition), 
  row.names = c("SRR6815098", "SRR6815093", "SRR6815094",  "SRR6815097", "SRR6815099", "SRR6815100", "SRR6815095", "SRR6815096")
)
#rownames(colData(ddsMat_SRP134109))

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Condition = c(IgG_control = "lightblue", Il1_ab = "orange")
)

# Make Heatmap with pheatmap function.
## See more in documentation for customization
pheatmap(mat = mat,  
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 6.5, # Make fonts smaller
         cellwidth = 20, # Make the cells wider
         show_colnames = F, cluster_cols = FALSE)

#Get fold change information
results_SRP134109_sig <- mapIds(x = org.Mm.eg.db,
                                            keys = rownames(results_SRP134109_sig),
                                            column = "SYMBOL",
                                            keytype = "ENSEMBL",
                                            multiVals = "first")

results_SRP134109_sig <- na.omit(results_SRP134109_sig)

genes_of_interest_listing_SRP134109 <- subset(results_SRP134109, rownames(results_SRP134109) %in% names_of_genes_of_interest_ensembl)
View(genes_of_interest_listing_SRP134109)


################################### Region data for SRP099054 project - Smooth muscle cell-specific deletion of Col15a1 #####################################

data_SRP099054 <- recount3::create_rse_manual(
  project = "SRP099054",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)

data_SRP099054

rowRanges(data_SRP099054)
colData(data_SRP099054)

## Find the project 
proj_info <- subset(
  mouse_projects,
  project == "SRP099054" & project_type == "data_sources"
)

## Create a RangedSummarizedExperiment (RSE) object at the gene level
rse_gene_SRP099054 <- create_rse(proj_info)

#save RSE object as rds file
saveRDS(rse_gene_SRP099054, file = "data_SRP099054/rse_gene_SRP099054.Rds")
rse_gene_SRP099054 <- readRDS("data_SRP099054/rse_gene_SRP099054.Rds")

assay(rse_gene_SRP099054, "counts") <- transform_counts(rse_gene_SRP099054)

countdata_SRP099054 <- read.delim("./data_SRP099054/GSE94661_collagen-atherosclerosis-in-vivo-rnaseq-rawcounts.txt", header = TRUE, row.names = 1)


# load Metadata information 
#Sample_title:	"in vivo_WT 1_SL142922"	"in vivo_WT 2_SL142923"	"in vivo_WT 3_SL142924"	"in vivo_KO 1_SL142925"	"in vivo_KO 2_SL142926"	"in vivo_KO 3_SL142927"
metadata_SRP099054 <- read.delim("./data_SRP099054/metadata.txt", row.names = 1)


#Deseq2 we compute log fold changes for the condition IgG vs. Il1b treatment
ddsMat_SRP099054 <- DESeqDataSetFromMatrix(countData = countdata_SRP099054,
                                           colData = metadata_SRP099054,
                                           design = ~condition)


# Find differential expressed genes
ddsMat_SRP099054 <- DESeq(ddsMat_SRP099054)

# Get results from testing with FDR adjust pvalues
results_SRP099054 <- results(ddsMat_SRP099054, pAdjustMethod = "fdr", alpha = 0.05)

#remove NA's
results_SRP099054 <- na.omit(results_SRP099054)

# Subset for only significant genes (q < 0.05)
results_SRP099054_sig <- subset(results_SRP099054, padj < 0.05)
results_SRP099054_sig <- na.omit(results_SRP099054_sig)
head(results_SRP099054_sig)

genes_of_interest_listing_SRP099054 <- subset(results_SRP099054, rownames(results_SRP099054) %in% names_of_genes_of_interest)
View(genes_of_interest_listing_SRP099054)
#Make heatmap showing these 5 genes of interest

#Heatmap
# Convert all samples to rlog
ddsMat_rlog_SRP099054 <- rlog(ddsMat_SRP099054, blind = FALSE)

mat <- assay(ddsMat_rlog_SRP099054[names_of_genes_of_interest])

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Condition = factor(colData(ddsMat_SRP099054)$condition), 
  row.names = c("SL142922", "SL142923", "SL142924",  "SL142925", "SL142926", "SL142927")
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Condition = c(KO = "lightblue", WT = "orange")
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

################################### Region data for SRP222892 project - Migratory and Dancing Macrophage Subsets in Atherosclerotic Lesions #####################################

data_SRP222892 <- recount3::create_rse_manual(
  project = "SRP222892",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)

data_SRP222892

rowRanges(data_SRP222892)
colData(data_SRP222892)

## Find the project 
proj_info <- subset(
  mouse_projects,
  project == "SRP222892" & project_type == "data_sources"
)

## Create a RangedSummarizedExperiment (RSE) object at the gene level
rse_gene_SRP222892 <- create_rse(proj_info)

#save RSE object as rds file
saveRDS(rse_gene_SRP222892, file = "data_SRP222892/rse_gene_SRP222892.Rds")
rse_gene_SRP222892 <- readRDS("data_SRP222892/rse_gene_SRP222892.Rds")

assay(rse_gene_SRP222892, "counts") <- transform_counts(rse_gene_SRP222892)

#Read countdata
countdata_SRP222892 <- as.matrix(fread("GSE137819_Counts.csv.gz"),rownames=1)

# load Metadata information 
#Sample_title:	"in vivo_WT 1_SL142922"	"in vivo_WT 2_SL142923"	"in vivo_WT 3_SL142924"	"in vivo_KO 1_SL142925"	"in vivo_KO 2_SL142926"	"in vivo_KO 3_SL142927"
metadata_SRP222892 <- read.delim("./data_SRP222892/metadata.txt", row.names = 1)


#Deseq2 we compute log fold changes for the condition IgG vs. Il1b treatment
ddsMat_SRP222892 <- DESeqDataSetFromMatrix(countData = countdata_SRP222892,
                                           colData = metadata_SRP222892,
                                           design = ~condition)


# Find differential expressed genes
ddsMat_SRP222892 <- DESeq(ddsMat_SRP222892)

# Get results from testing with FDR adjust pvalues
results_SRP222892 <- results(ddsMat_SRP222892, pAdjustMethod = "fdr", alpha = 0.05)

#remove NA's
results_SRP222892 <- na.omit(results_SRP222892)

# Subset for only significant genes (q < 0.05)
results_SRP222892_sig <- subset(results_SRP222892, padj < 0.05)
results_SRP222892_sig <- na.omit(results_SRP222892_sig)
head(results_SRP222892_sig)

genes_of_interest_listing_SRP222892 <- subset(results_SRP222892, rownames(results_SRP222892) %in% names_of_genes_of_interest)
View(genes_of_interest_listing_SRP222892)

#Make heatmap showing these 5 genes of interest
# Convert all samples to rlog
ddsMat_rlog_SRP222892 <- rlog(ddsMat_SRP222892, blind = FALSE)

mat <- assay(ddsMat_rlog_SRP222892[names_of_genes_of_interest])

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Condition = factor(colData(ddsMat_SRP222892)$condition), 
  Replicate = factor(colData(ddsMat_rlog_SRP222892)$replicate),
  row.names = c("GFP_1", "GFP_2", "GFP_3", "GFP_4", "GFP_5", "GFP_6", "DP_1" , "DP_2",  "DP_4",  "DP_5",  "YFP_1", "YFP_3", "YFP_4", "YFP_5", "YFP_6", "DN_1",  "DN_3",  "DN_6")
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Condition = c(DN = "lightblue", DP = "orange", GFP = "green", YFP = "grey")
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



#########################################  Region data for SRP119232 project - Mus musculus Transcriptome or Gene expression ##################

data_SRP119232 <- recount3::create_rse_manual(
  project = "SRP119232",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)

data_SRP119232

rowRanges(data_SRP119232)
colData(data_SRP119232)

## Find the project 
proj_info <- subset(
  mouse_projects,
  project == "SRP119232" & project_type == "data_sources"
)

## Create a RangedSummarizedExperiment (RSE) object at the gene level
rse_gene_SRP119232 <- create_rse(proj_info)

#save RSE object as rds file
saveRDS(rse_gene_SRP119232, file = "data_SRP119232/rse_gene_SRP119232.Rds")
rse_gene_SRP119232 <- readRDS("data_SRP119232/rse_gene_SRP099054.Rds")

assay(rse_gene_SRP119232, "counts") <- transform_counts(rse_gene_SRP119232)

#Write sra.sample_description to file
write.table(rse_gene_SRP119232$sra.sample_description, "./data_SRP119232/sample_description.tsv")

#abort -> dataset is far too big
