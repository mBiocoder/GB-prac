library(data.table)
library(ggplot2)
library(magrittr)
library(dplyr)

files <- list.files('extdata', full.names = T)
names(files) <- basename(files)
tables <- lapply(files, fread)
dt <- rbindlist(tables, idcol = 'filepath')


ggplot(dt, aes(x=max_skipped_exon, fill=filepath)) + geom_bar() + scale_y_log10()
ggplot(dt, aes(x=max_skipped_exon)) + geom_bar()
ggplot(dt, aes(x=max_skipped_exon)) + geom_bar() + scale_y_log10()
ggplot(dt, aes(x=max_skipped_exon)) + geom_histogram() + scale_y_log10() + facet_wrap(vars(filepath))

ggplot(dt, aes(x=max_skipped_exon, y=filepath)) + geom_boxplot()
ggplot(dt, aes(x=max_skipped_exon, color=filepath)) + geom_histogram()
ggplot(dt, aes(x=max_skipped_exon, color=filepath)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") + scale_y_log10()

# cumulative plot exons!!!
ggplot(dt %>% 
         group_by(filepath) %>% 
         arrange(max_skipped_exon) %>% 
         mutate(rn = row_number())) + 
  geom_step(aes(x=max_skipped_exon, y=rn, color=filepath)) + labs(x="number of max skipped exons", y="cumulative number of appearance")






ggplot(dt, aes(x=max_skipped_bases, fill=filepath)) + geom_bar() + scale_y_log10()
ggplot(dt, aes(x=max_skipped_bases)) + geom_bar() + scale_y_log10()

ggplot(dt, aes(x=max_skipped_bases, color=filepath)) + geom_histogram()
ggplot(dt, aes(x=max_skipped_bases, color=filepath))+geom_histogram(aes(y=cumsum(..count..)))
ggplot(dt, aes(x=max_skipped_bases, color=filepath))+ stat_ecdf()

# cumulative plot bases!!!
ggplot(dt %>% 
         group_by(filepath) %>% 
         arrange(max_skipped_bases) %>% 
         mutate(rn = row_number())) + 
  geom_step(aes(x=max_skipped_bases, y=rn, color=filepath))

ggplot(dt, aes(x=max_skipped_bases, y=filepath)) + geom_boxplot()

ggplot(dt, aes(x=max_skipped_bases)) + geom_bar() + scale_y_log10() + facet_wrap(vars(filepath))

