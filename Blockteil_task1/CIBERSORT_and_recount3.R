#BiocManager::install(c("DESeq2", "edgeR", "RUVSeq", "recount3"))

## downloading gtex data with recount3
library(recount3)
library(DESeq2)
library(edgeR)
library(RUVSeq)

## get table of available projects from recount3
human_projects <- available_projects()

dim(human_projects)
head(human_projects)

ds <- recount3::create_rse_manual(project = "SRP082563",
                                  project_home = "data_sources/sra",
                                  organism = "mouse",
                                  annotation = "gencode_v23",
                                  type = "gene")

head(ds)

## Find all the mouse projects
mouse_projects <- available_projects(organism = "mouse")

## Explore the results
dim(mouse_projects)
head(mouse_projects)



###############################
## find gtex projects
gtex <- subset(human_projects, file_source == "gtex" & project_type == "data_sources")
gtex

## gtex data in recount3 is stored according to tissue type
## we want the blood vessel data in this case
proj_info <- subset(gtex, project == "BLOOD_VESSEL" & project_type == "data_sources")

## make the summarized experiment object from the blood vessel gtex data
rse_blood_vessel <- create_rse(proj_info)

## we only want the coronary artery samples, stored in gtex.smtsd as "Artery - Coronary"
rse_coronary <- rse_blood_vessel[ , rse_blood_vessel$gtex.smtsd == "Artery - Coronary"]

## transform counts as with the reference data
assay(rse_coronary, "counts") <- transform_counts(rse_coronary)
saveRDS(rse_coronary, file = "CIBERSORT/rse_coronary.Rds")



#################################### Preprocessing ############################################

## read in reference data for batch correction
rse <- readRDS(file = "CIBERSORT/rse_coronary.Rds")

## filter out low counts to ease computation
low_counts <- apply(assay(rse, "counts"), MARGIN = 1, FUN = function(x){all(x < 1000)})
rse <- rse[-which(low_counts), ]

## create expression set for edgeR
celltypes <- rse$celltype
set <- newSeqExpressionSet(assay(rse, "counts"), 
                           phenoData = data.frame(celltypes, row.names = colnames(rse)))
design <- model.matrix(~celltypes, data = pData(set))
genes <- rownames(rse)

## run edgeR
y <- DGEList(counts = counts(set), group = celltypes) %>% 
  calcNormFactors(object = ., method = "TMM") %>%
  estimateGLMCommonDisp(y = ., design = design) %>%
  estimateGLMTagwiseDisp(y = ., design = design)
fit <- glmFit(y = y, design = design)
res <- residuals(fit, type = "pearson")

## run RUV seq with RUVr method for blind batch correction
set0 <- RUVr(set, genes, k = 20, res)

## create summarized experiment object for cibersort preparation
## see ciber_prep.R for further steps
ruv <- SummarizedExperiment(assays = list(counts = set0@assayData$normalizedCounts), 
                            rowData = rowData(rse), colData = colData(rse))
saveRDS(ruv, file = "data/ruv_reference.Rds")