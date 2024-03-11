x <- "probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy"

df <- data.table::fread(x)
# View(df)

finaldf <- data.frame(matrix(ncol = 6, nrow = 0))


for(j in 1:nrow(df)) {

  probe <- unlist(df[j,1])
  rest <- unlist(df[j,3:6])
  
  splitted <- strsplit(unlist(df[j,2]),",")[[1]]
  
  for(i in 1:length(splitted)) {
    newrow <- c(probe, splitted[i], rest)
    finaldf <- rbind(finaldf, newrow)
  }

}

colnames(finaldf) <- colnames(df)

finalname <- paste(substring(x, 1, (nchar(x) - 4)), "_Modified.csv", sep = "")
write.table(finaldf, finalname, sep = ",", row.names = FALSE, col.names = TRUE)