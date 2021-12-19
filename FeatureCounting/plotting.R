library(data.table)
library(ggplot2)

ebna <- fread("../plottingData/ebna_hisat.plotting")
hes <- fread("../plottingData/hes_star.plotting")
nookaew <- fread("../plottingData/nookaew_cm.plotting")


#### nookaew ###########
ggplot(nookaew, aes(x=log2rpkm)) + geom_histogram(bins=100, aes(y=..density..), color="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  labs(x="log2 FPKM value", y="density", title="FPKM value distribution (all reads)", subtitle = "for nookaew_cm")

ggplot(nookaew, aes(x=log2rpkm_pcrZero)) + geom_histogram(bins=100, aes(y=..density..), color="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  labs(x="log2 FPKM value", y="density", title="FPKM value distribution (only non-PCR-duplicated reads)", subtitle = "for nookaew_cm")

ggplot(nookaew) + geom_density(aes(x=log2rpkm, fill="red", color="red"), alpha=.3) +
  geom_density(aes(x=log2rpkm_pcrZero, fill="blue", color="blue"), alpha = .3)+guides(colour=FALSE)+
  labs(x="log2 FPKM value", y="density", title="FPKM value distribution", subtitle = "for nookaew_cm") + 
  scale_fill_discrete(name="FPKM of", labels=c("all reads", "non-PCR-duplicated reads"))

#### hes ###########
ggplot(hes, aes(x=log2rpkm)) + geom_histogram(bins=100, aes(y=..density..), color="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  labs(x="log2 FPKM value", y="density", title="FPKM value distribution (all reads)", subtitle = "for hes_star")

ggplot(hes, aes(x=log2rpkm_pcrZero)) + geom_histogram(bins=100, aes(y=..density..), color="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  labs(x="log2 FPKM value", y="density", title="FPKM value distribution (only non-PCR-duplicated reads)", subtitle = "for hes_star")

ggplot(hes) + geom_density(aes(x=log2rpkm, fill="red", color="red"), alpha=.3) +
  geom_density(aes(x=log2rpkm_pcrZero, fill="blue", color="blue"), alpha = .3)+guides(colour=FALSE)+
  labs(x="log2 FPKM value", y="density", title="FPKM value distribution", subtitle = "for hes_star") + 
  scale_fill_discrete(name="FPKM of", labels=c("all reads", "non-PCR-duplicated reads"))

#### ebna ###########
ggplot(ebna, aes(x=log2rpkm)) + geom_histogram(bins=100, aes(y=..density..), color="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  labs(x="log2 FPKM value", y="density", title="FPKM value distribution (all reads)", subtitle = "for ebna_hisat")

ggplot(ebna, aes(x=log2rpkm_pcrZero)) + geom_histogram(bins=100, aes(y=..density..), color="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  labs(x="log2 FPKM value", y="density", title="FPKM value distribution (only non-PCR-duplicated reads)", subtitle = "for ebna_hisat")

ggplot(ebna) + geom_density(aes(x=log2rpkm, fill="red", color="red"), alpha=.3) +
  geom_density(aes(x=log2rpkm_pcrZero, fill="blue", color="blue"), alpha = .3)+guides(colour=FALSE)+
  labs(x="log2 FPKM value", y="density", title="FPKM value distribution", subtitle = "for ebna_hisat") + 
  scale_fill_discrete(name="FPKM of", labels=c("all reads", "non-PCR-duplicated reads"))

#### all ############
ebna[, file:="ebna"]
hes[, file:="hes"]
nookaew[, file:="nookaew"]
dt <- rbind(ebna, hes, nookaew)

ggplot(dt, aes(x=log2rpkm, fill=file)) + geom_density(alpha=.3) +
  labs(x="log2 FPKM value", y="density", title="FPKM value distribution (all reads)", subtitle = "comparison of all three bam files") + 
  scale_fill_discrete(name="bam file", labels=c("ebna_hisat", "hes_star", "nookaew_cm"))

#ggplot(dt) + geom_density(aes(x=rpkm, fill="red", color="red"), alpha=.3) + scale_x_log10() +
#  geom_density(aes(x=rpkm_pcrZero, fill="blue", color="blue"), alpha = .3)+guides(colour=FALSE)+
#  labs(x="log10 FPKM value", y="density", title="FPKM value distribution", subtitle = "for nookaew_cm") + 
#  scale_fill_discrete(name="FPKM of", labels=c("all reads", "non-PCR-duplicated reads"))


#ggplot(dt) + geom_histogram(aes(x=rpkm, fill="red"), bins=200, alpha=.5, position="identity")+
#  scale_x_log10() + geom_histogram(aes(x=rpkm_pcrZero, fill="blue"), bins=200, alpha=.5, position="identity")
