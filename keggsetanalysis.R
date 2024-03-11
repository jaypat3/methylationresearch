Hypergeometric <- function(data1,data2,n) {
  t <- length(intersect(data1,data2))
  if(length(data1) >= length(data2)) {
    a <- length(data1)
    b <- length(data2)
  }
  else {
    a <- length(data2)
    b <- length(data1)
  }
  toreturn <- vector()
  toreturn <- c(t,sum(dhyper(t:b, a, n-a, b)))
  return(toreturn)
}

datasetoptions <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC")
perc <- "10Perc.csv"

for(whichdata in datasetoptions) {
  datasetname <- paste(paste(paste(paste(whichdata,"datasets/GeneTable",sep=""),whichdata,sep=""),"PValsMOD_",sep=""),perc,sep="")
  datasettype <- gsub(".*Table(.*)\\PVals.*", "\\1", datasetname)
  
  
  df <- data.table::fread(datasetname)
  
  
  tempdf <- data.frame(matrix(ncol=4,nrow=0))
  
  
  
  rownamevec <- vector()
  
  for(filename in list.files("kegggenesets/")) {
    location <- paste("kegggenesets/",filename,sep="")
    genesetname <- gsub(".*KEGG_(.*)\\.v2023.*","\\1", location)
    set1 <- unlist(data.table::fread(location))
    rownamevec <- c(rownamevec,genesetname)
    
    gtmtset <- unlist(unique(df[df$GTMT <= 0.05 & df$GTMTnSig > 0.05 & df$GTCG > 0.05 & df$GTCD > 0.05 & df$GTMut > 0.05]$Gene))
    gtmbset <- unlist(unique(df[df$GTMB <= 0.05 & df$GTMBnSig > 0.05 & df$GTCG > 0.05 & df$GTCD > 0.05 & df$GTMut > 0.05]$Gene))
    gbmtset <- unlist(unique(df[df$GBMT <= 0.05 & df$GBMTnSig > 0.05 & df$GBCG > 0.05 & df$GBCD > 0.05 & df$GBMut > 0.05]$Gene))
    gbmbset <- unlist(unique(df[df$GBMB <= 0.05 & df$GBMBnSig > 0.05 & df$GBCG > 0.05 & df$GBCD > 0.05 & df$GBMut > 0.05]$Gene))
    
    
    gtmtoutput <- Hypergeometric(set1,gtmtset,nrow(df))[[2]]
    gtmboutput <- Hypergeometric(set1,gtmbset,nrow(df))[[2]]
    gbmtoutput <- Hypergeometric(set1,gbmtset,nrow(df))[[2]]
    gbmboutput <- Hypergeometric(set1,gbmbset,nrow(df))[[2]]
    
    setvec <- c(gtmtoutput,gtmboutput,gbmtoutput,gbmboutput)
    tempdf <- rbind(tempdf,setvec)
  }
  
  pvalprechange <- unlist(tempdf)
  pvalchange <- p.adjust(pvalprechange,method="fdr")
  finaldf <- data.frame(matrix(pvalchange,ncol=4,nrow=nrow(tempdf)))
  colnames(finaldf) = c("GTMT", "GTMB", "GBMT", "GBMB")
  rownames(finaldf) = rownamevec
  
  write.table(finaldf, paste(datasettype, "_10perc_genesetenrichment.csv",sep=""),sep=",", row.names=TRUE, col.names = TRUE)
}