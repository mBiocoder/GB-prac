#!/usr/bin/Rscript --vanilla

#files <- c("sample1.psi", "sample2.psi","sample3.psi","sample4.psi","sample5.psi","sample6.psi","sample7.psi","sample8.psi","sample9.psi", "sample10.psi")
#group <- c(1,1,1,1,1,2,2,2,2,2)
#r <- diff.splicing(files, group)
#write.table(r, file="diff_analysis_results.tsv", sep='\t', quote=FALSE, row.names = FALSE)

#i <- c(5, 14, 13, 15, 29, 8, 12, 13, 22, 43)
#names(i) <- c("sample1", "sample2","sample3","sample4","sample5","sample6","sample7","sample8","sample9","sample10")
#t <- c(51, 61, 52, 49, 51, 58, 52, 43, 57, 64)
#names(t) <- c("sample1", "sample2","sample3","sample4","sample5","sample6","sample7","sample8","sample9","sample10")
#g <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
#LRS(i,t,g)


LRS <- function(incl, total, group){
  # change the function so that it only takes single rows as input...
  p0 <- sum(incl)/sum(total)
  p1 <- sum(incl[group == 1])/sum(total[group == 1])
  p2 <- sum(incl[group == 2])/sum(total[group == 2])
  
  # compute likelihoods:
  lreduced <- prod(choose(total, incl)*(p0^incl)*((1-p0)^(total-incl)))
    
  full <- prod(choose(total[group == 1], incl[group == 1])*
                   (p1^incl[group == 1])*
                   ((1-p1)^(total[group == 1]-incl[group == 1])))
    
  lfull <- full * prod(choose(total[group == 2], incl[group == 2])*
                           (p2^incl[group == 2])*
                           ((1-p2)^(total[group == 2]-incl[group == 2])))
  
  
  lrs <- (-2)*log(lreduced/lfull)
  
  pvalue <- pchisq(lrs, df = 1,lower.tail = FALSE)
  
  return(list(p0=p0, p1=p1, p2=p2, llreduced = log(lreduced), llfull = log(lfull), lrs = lrs, pvalue = pvalue))
}


diff.splicing <- function(psi.files, group){
  result <- read.delim(psi.files[1])[1:2]
  row <- c()
  for(x in 1:nrow(result)){
    row <- c(row, paste0("exon", x))
  }
  
  result <- cbind(result, row)
  inclusion <- data.frame(result$row)
  colnames(inclusion) <- c("row")
  total <- data.frame(result$row)
  colnames(total) <- c("row")
  
  for(x in 1:length(psi.files)){
    file <- psi.files[x]
    data <- read.delim(file)
    data <- merge(data, result, by=c("gene", "exon"), all = TRUE)
    colnames(data)[3] <- paste0("sample",x)
    colnames(data)[5] <- paste0("sample",x)
    # remove possible NA because we used merge with all
    data[is.na(data)] <- 0
    inclusion <- merge(inclusion, data[c(3,7)], by=c("row"))
    total <- merge(total, data[c(5,7)], by=c("row"))
  }
  
  row_names <- as.vector(inclusion[,1])
  inclusion <- data.matrix(inclusion[2:ncol(inclusion)], rownames.force = TRUE)
  rownames(inclusion) <- row_names
  
  row_names <- as.vector(total[,1])
  total <- data.matrix(total[2:ncol(total)], rownames.force = TRUE)
  rownames(total) <- row_names
  
  lrs_result <- data.frame(LRS(inclusion[1,], total[1,], group))
  lrs_result <- cbind(lrs_result, c(rownames(total)[1]))
  colnames(lrs_result)[8] <- "row"
  
  
  for(i in 2:nrow(total)){
    lrs_frame <- data.frame(LRS(inclusion[i,], total[i,], group))
    lrs_frame <- cbind(lrs_frame, rownames(total)[i])
    colnames(lrs_frame)[8] <- "row"
    lrs_result <- rbind(lrs_result, lrs_frame)
  }
  
  lrs_result <- merge(lrs_result, result, by = c("row"))
  lrs_result <- lrs_result[, -which(names(lrs_result) %in% c("row"))]
  lrs_result <- lrs_result[, c(8, 9, 1, 2, 3, 4, 5, 6, 7)]
  
  padj <- p.adjust(lrs_result$pvalue, method = "BH")
  lrs_result <- cbind(lrs_result, padj)
  
  return(lrs_result)
}