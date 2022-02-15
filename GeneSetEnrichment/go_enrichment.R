library(ggplot2)
library(data.table)
library(dplyr)
library(plyr)

#################################### Read files in ################################################

#input files
goa_human_ensembl <- as.data.frame(fread("./goa_human_ensembl.tsv"))
simul_exp_go_bp_ensembl <- as.data.frame(fread("./simul_exp_go_bp_ensembl.tsv"))

#result files between minsize and maxsize
overlap_out <- read.table("./overlapout.tsv", sep = '\t', header = TRUE)
out <- as.data.frame(fread("./out.tsv"))
#out <- subset(out, select = -c(shortest_path_to_a_true)) #nötig, da extra spalte nach dem einlesen entsteht
#colnames(out) <- c("term", "name", "size", "is_true", "noverlap", "hg_pval", "hg_fdr", "fej_pval", "fej_fdr", "ks_stat", "ks_pval", "ks_fdr", "count_signif", "shortest_path_to_a_true")

#result files for all genesets
overlap_out_all <- read.table("./overlapout_all.tsv", sep = '\t', header = TRUE)
out_all <- as.data.frame(fread("./out_all.tsv"))

#SGs between minsize and maxsize in enrichmentout
sgInGenesets <- overlap_out_all <- read.table("./sgInGenesets.tsv", sep = '\t', header = TRUE)

#SGs between minsize and maxsize in ovelapout
sgInOverlapGenesets <- overlap_out_all <- read.table("./sgInOverlapGenesets.tsv", sep = '\t', header = TRUE)

################################################################################################

#list of genes which are significant 
all_signif_genes <- simul_exp_go_bp_ensembl %>% filter(simul_exp_go_bp_ensembl$signif == "TRUE")
write.table(all_signif_genes, file='all_significant_genes.tsv', quote=FALSE, sep='\t')

# list of simulated_enriched_GO_ids
simulated_enriched_GO_ids <- c("GO:0098743","GO:0098744","GO:0098745","GO:0098746",
                               "GO:0098747","GO:0098748","GO:0098749","GO:0098750", 
                               "GO:0098751","GO:0098752","GO:0098753","GO:0098754",
                               "GO:0098755","GO:0098756","GO:0098757","GO:0098758")

######################################### Plotting ######################################################

#scatter plot of p-value against size
ggplot(out, aes(x=size, y=hg_pval)) + geom_point() 
ggplot(out, aes(x=size, y=fej_pval)) + geom_point() #use this for comparison with DAVID
ggplot(out, aes(x=size, y=ks_pval)) + geom_point() 

#scatterplot of num_overlap against path_length
ggplot(overlap_out, aes(x=path_length, y=num_overlapping)) + geom_point() 

#scatterplot noverlap against size from out file
p <- ggplot(out, aes(x=size, y=noverlap)) + geom_point()
p +  geom_smooth(method='lm', formula= y~x)

# create multiple linear model
lm_fit <- lm(noverlap ~ size, data=out)
summary(lm_fit)

library("ggpubr")
ggscatter(out, x = "size", y = "noverlap", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson")

#scatterplot score against size
ggplot(out, aes(x=ks_stat, y=size)) + geom_point()

#p-value histograms to show distribution of p-values
ggplot(out, aes(x=hg_pval)) + geom_histogram()
ggplot(out, aes(x=fej_pval)) + geom_histogram()
ggplot(out, aes(x=ks_pval)) + geom_histogram()

#p-value densities distribution of p-values
ggplot(out, aes(x=hg_pval)) + geom_density()
ggplot(out, aes(x=fej_pval)) + geom_density()
ggplot(out, aes(x=ks_pval)) + geom_density()

#ES histogram to show distribution of enrichment scores
ggplot(out, aes(x=ks_stat)) + geom_histogram()

#ES densities to show distribution of enrichment scores
ggplot(out, aes(x=ks_stat)) + geom_density() 


#Plot distribution and histogram of significant genes for gene sets between minsize and maxsize: mean= 43.79228
ggplot(out, aes(x=count_signif)) + geom_density() + labs(x="Significant genes") + geom_vline(aes(xintercept=mean(count_signif)),
                                                                                     color="blue", linetype="dashed", size=1)
ggplot(out, aes(x=count_signif)) + geom_histogram() + labs(x="Significant genes")


#Plot distribution and histogram of significant genes for all gene sets: mean= 8.792615
ggplot(out, aes(x=count_signif)) + geom_density() + labs(x="Significant genes") + geom_vline(aes(xintercept=mean(count_signif)),
                                                                                             color="blue", linetype="dashed", size=1)
ggplot(out_all, aes(x=count_signif)) + geom_histogram() + labs(x="Significant genes") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),limits=c(0, 100000), breaks=c(0, 1, 10, 100, 1000, 10000))


#Inside minmax size
#Sorted tables of gene sets according to number of SGs, score and p-value; high number of SGs and score and low p-value

table_minmax_SG <- arrange(out,desc(out$count_signif))
View(table_minmax_SG)

table_minmax_score <- arrange(out,desc(out$ks_stat))
View(table_minmax_score)

table_minmax_fej_pval <- arrange(out,out$fej_pval)
View(table_minmax_fej_pval)

table_minmax_hg_pval <- arrange(out,out$hg_pval)
View(table_minmax_fej_pval)

table_minmax_ks_pval <- arrange(out,out$ks_pval)
View(table_minmax_ks_pval)


#All genesets
table_all_SG <- arrange(out_all,desc(out_all$count_signif))
View(table_all_SG)

table_all_score <- arrange(out_all,desc(out_all$ks_stat))
View(table_all_score)

table_all_fej_pval <- arrange(out_all,out_all$fej_pval)
View(table_all_fej_pval)

table_all_hg_pval <- arrange(out_all,out_all$hg_pval)
View(table_all_hg_pval)

table_all_ks_pval <- arrange(out_all,out_all$ks_pval)
View(table_all_ks_pval)


#Plot distirubution of SGs in how many gene sets
ggplot(sgInGenesets, aes(x=count)) + geom_density() + labs(x="Significant genes in gene sets")
#Are there unique SGs in only one set? -> Answer: TRUE, e.g. TEF, ASTN2, RSC1A1, etc.
#any(sgInGenesets$count == 1)

#Plot distirubution of SGs in how many overlapping gene sets
ggplot(sgInOverlapGenesets, aes(x=count)) + geom_density() + labs(x="Significant genes in overlapping gene sets")
#Are there unique SGs in only one set? -> Answer: TRUE, e.g. GMIP, CIPC, PASD1, etc.
any(sgInOverlapGenesets$count == 1) 

