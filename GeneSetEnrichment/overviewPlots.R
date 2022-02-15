library(ggplot2)
library(data.table)
library(dplyr)
library(qpcR)

dt <- fread("biological_process.leavesToRoot.txt")
colnames(dt)[2] <- "biological_process"
dt <- merge(dt, fread("cellular_component.leavesToRoot.txt"), by="set", all = TRUE)
colnames(dt)[3] <- "cellular_component"
dt <- merge(dt, fread("molecular_function.leavesToRoot.txt"), by="set", all = TRUE)
colnames(dt)[4] <- "molecular_function"
dt <- dt[, c(2, 3, 4)]
dt <- melt(dt, na.rm=TRUE, value.name="leaves_to_root", variable.name = "namespace")


ggplot(dt, aes(x = leaves_to_root, fill = namespace)) +
  geom_histogram(position = "identity", alpha = 0.4, binwidth = 1) +
  scale_x_continuous(breaks=0:12) +
  labs(title="Distribution of shortest paths from leaf to root", x="path length")


dt2 <- fread("biological_process.numGenesInSets.txt")
colnames(dt2)[2] <- "biological_process"
dt2 <- merge(dt2, fread("cellular_component.numGenesInSets.txt"), by="set", all = TRUE)
colnames(dt2)[3] <- "cellular_component"
dt2 <- merge(dt2, fread("molecular_function.numGenesInSets.txt"), by="set", all = TRUE)
colnames(dt2)[4] <- "molecular_function"
dt2 <- dt2[, c(2, 3, 4)]
dt2 <- melt(dt2, na.rm=TRUE, value.name="set_sizes", variable.name = "namespace")

ggplot(dt2, aes(x = set_sizes, color = namespace)) + stat_ecdf() +
  scale_y_log10()

ggplot(dt2, aes(x = set_sizes, fill = namespace)) +
  geom_histogram(position = "identity", alpha = 0.4)

# only values between min and maxsize
ggplot(dt2[set_sizes >= 0 & set_sizes <= 10,], aes(x = set_sizes, fill = namespace)) +
  geom_histogram(position = "identity", alpha = 0.4, binwidth=1) +
  scale_x_continuous(breaks=0:10)+
  labs(title="Distribution of gene set sizes between 0 and 10", x="size")

ggplot(dt2[set_sizes >= 10 & set_sizes <= 50,], aes(x = set_sizes, fill = namespace)) +
  geom_histogram(position = "identity", alpha = 0.4, binwidth=1)+
  labs(title="Distribution of gene set sizes between 10 and 50", x="size")

ggplot(dt2[set_sizes >= 50 & set_sizes <= 500,], aes(x = set_sizes, fill = namespace)) +
  geom_histogram(position = "identity", alpha = 0.4) +
  labs(title="Distribution of gene set sizes between 50 and 500", x="size")

ggplot(dt2[set_sizes >= 500 & set_sizes <= 5000, ], aes(x = set_sizes, fill = namespace)) +
  geom_histogram(position = "identity", alpha = 0.4, binwidth = 200)

ggplot(dt2[set_sizes >= 5000, ], aes(x = set_sizes, fill = namespace)) +
  geom_histogram(position = "identity", alpha = 0.4, binwidth = 200)





dt3 <- fread("biological_process.sizeDiffs.txt")
colnames(dt3)[1] <- "biological_process"
dt3 <- qpcR:::cbind.na(dt3, fread("cellular_component.sizeDiffs.txt"))
colnames(dt3)[2] <- "cellular_component"
dt3 <- qpcR:::cbind.na(dt3, fread("molecular_function.sizeDiffs.txt"))
colnames(dt3)[3] <- "molecular_function"
dt3 <- melt(dt3, na.rm=TRUE, value.name="size_diffs", variable.name = "namespace")

ggplot(dt3, aes(x = size_diffs, fill = namespace)) +
  geom_histogram(position = "identity", alpha = 0.4, binwidth = 500) +
  scale_y_log10() + labs(title="Distribution of size differences between children and parents",
                         x="set size difference")


ggplot(dt3[order(dt3[,size_diffs]),],
       aes(x = size_diffs, y = ave(namespace == namespace, namespace, FUN = cumsum),
           col = namespace)) + geom_step() +
  labs(title="Size differences between parents and children",
       x="set size difference", y="cumulative count")

