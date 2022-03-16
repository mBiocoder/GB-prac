library(maSigPro)
library(data.table)
library(ggplot2)

################################## Multi-time series data #################################

data(data.abiotic)
data(edesign.abiotic)
edesign.abiotic

head(data.abiotic)
colnames(data.abiotic)
rownames(edesign.abiotic)
colnames(edesign.abiotic)
rownames(data.abiotic)


#defining the regression model
design <- make.design.matrix(edesign.abiotic, degree = 2)
design$groups.vector

#finding significant genes
fit <- p.vector(data.abiotic, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)

fit$i # returns the number of significant genes
fit$alfa # gives p-value at the Q false discovery control level
fit$SELEC # is a matrix with the significant genes and their expression values


#finding significant differences
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)

#obtain list of significant genes
sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
names(sigs)


names(sigs$sig.genes)
names(sigs$sig.genes$ColdvsControl)


#Venn diagrams
suma2Venn(sigs$summary[, c(2:4)])
suma2Venn(sigs$summary[, c(1:4)])

#see specific genes from venn diagram comparisons
sigs$sig.genes$SaltvsControl$g

see.genes(sigs$sig.genes$ColdvsControl, show.fit = T, dis =design$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 9)


#Plot groups
STMDE66 <- data.abiotic[rownames(data.abiotic)=="STMDE66", ]
PlotGroups (STMDE66, edesign = edesign.abiotic)

#add regression curve to plot
PlotGroups (STMDE66, edesign = edesign.abiotic, show.fit = T,
            + dis = design$dis, groups.vector = design$groups.vector)


############################################## Dataset 5 ####################################

#import countdata with normalized counts
countdata <- read.table("./CIBERSORT/D5_mix.tsv", header = TRUE, row.names = 1)
countdata <- na.omit(countdata)
head(countdata)

# Import metadata file
# - make row names the matching sampleID's from the countdata
metadata <- read.delim("./sample.txt", row.names = 1)

#remove taxid and condition column
metadata$taxid <- NULL
metadata$condition <- NULL

#Add time column to mapping file
metadata$Time <- c(0,0,0,0,1,1,1,1,2,2,2,2,4,4,4,4)
#metadata$Time <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)

#Add replicate column to mapping file
metadata$Replicate <- c(1, 2, 3, 4, 1, 2, 3, 4,1, 2, 3, 4,1, 2, 3, 4)

#Add control column to mapping file
metadata$Control <- c(1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)

#Add timepoint 1 as condition column to mapping file
metadata$Timepoint_1 <- c(0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0)

#Add timepoint 2 as condition column to mapping file
metadata$Timepoint_2 <- c(0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0)

#Add timepoint 4 as condition column to mapping file
metadata$Timepoint_4 <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1)

# Make sure ID's are correct
metadata
colnames(metadata)
rownames(metadata)


##defining the regression model
design <- make.design.matrix(metadata, degree = 2)
design$groups.vector

#finding significant genes
fit <- p.vector(countdata, design, Q = 0.05, MT.adjust = "BH")

fit$i # returns the number of significant genes: 6199 -> 1594 -> 7054
fit$alfa # gives p-value at the Q false discovery control level
fit$SELEC # is a matrix with the significant genes and their expression values


#finding significant differences
tstep <- T.fit(fit, step.method = "forward", alfa = 0.05)

get<-get.siggenes(tstep, vars="groups")
names(get)
names(get$sig.genes)
names(get$sig.genes$Timepoint_1vsControl)
get$summary

#Venn diagrams
suma2Venn(get$summary[, c(1:4)])
suma2Venn(get$summary[, c(2:4)])
#suma2Venn(get$summary[, c(3,4)])


#see specific genes from venn diagram comparisons
get$sig.genes$Timepoint_4vsControl$g

see.genes(get$sig.genes$Timepoint_4vsControl)
see.genes(get$sig.genes$Timepoint_4vsControl, show.fit = T, dis =design$dis)

see.genes(get$sig.genes$Timepoint_1vsControl, show.fit = T, dis =design$dis,
          cluster.method="hclust", k = 9)




####################################################################################
#Plot groups
gene_table <- countdata[rownames(countdata)=="Tnfaip8", ]
gene_table

#melt table from wide to long
df1 <- melt(gene_table, id.vars = 0, measure.vars = c("WT_0h_rep_1", "WT_0h_rep_2", "WT_0h_rep_3", "WT_0h_rep_4",
                                               "WT_1h_rep_1", "WT_1h_rep_2", "WT_1h_rep_3", "WT_1h_rep_4",
                                               "WT_2h_rep_1", "WT_2h_rep_2", "WT_2h_rep_3", "WT_2h_rep_4",
                                               "WT_4h_rep_1", "WT_4h_rep_2", "WT_4h_rep_3", "WT_4h_rep_4"),
     variable.name = "sample", value.name = "counts") 

#add time column to df1
df1$time <- c(0,0,0,0,1,1,1,1,2,2,2,2,4,4,4,4)
df1$time <- as.character(df1$time)
#Add replicate column to df1
df1$replicate <- c(1, 2, 3, 4, 1, 2, 3, 4,1, 2, 3, 4,1, 2, 3, 4)
df1$replicate <- as.character(df1$replicate)

summary(df1)

#Plot for gene time progression
ggplot(df1, aes(x=factor(time), y=counts, group=replicate)) +
  geom_line(aes(color=replicate))+
  geom_point(aes(color=replicate)) + xlab("time") + 
  ggtitle("Tnfaip8") + theme(plot.title = element_text(hjust = 0.5)) 

