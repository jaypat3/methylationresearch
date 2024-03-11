#This will construct Gene tables and Atypical Probe Tables for TCGA Datasets

DATASET <- "BRCA"

LoadLibraries <- function() {
  library(FedData)
  library(r2r)
}

SeparateTumorsAndNormals <- function(df1) {
  x <- colnames(df1)
  y <- substr_right(x, 3)
  z <- which(y == "-11")
  
  returnList <- list()
  returnList[[1]] <- df1[,!..z]
  returnList[[2]] <- df1[,..z]
  return(returnList)
}

FindLocation <- function(probeStart, chromStart, chromEnd) {
  if(length(chromStart) == 0) {
    return("NA")
  }
  if(chromStart <= probeStart && probeStart < chromEnd) return("GeneBody")
  else if(chromEnd <= probeStart) return("AfterGene")
  else if(probeStart < chromStart && (chromStart - probeStart) <= 1000) return("Promoter")
  else return("BeforeGene")
}

PrepareGeneSamples <- function(GeneExpTumordf, mainGene) {
  GeneExpRow <- GeneExpTumordf[which(GeneExpTumordf[,1] == mainGene),-1]
  GeneTop20 <- quantile(GeneExpRow, probs = c(0.8),na.rm = TRUE)
  GeneBot20 <- quantile(GeneExpRow, probs = c(0.2),na.rm = TRUE)
  GeneExpTop <- colnames(GeneExpTumordf)[which(GeneExpRow > GeneTop20)+1]
  GeneExpBottom <- colnames(GeneExpTumordf)[which(GeneExpRow < GeneBot20)+1]
  toExport <- list()
  toExport[[1]] <- GeneExpTop
  toExport[[2]] <- GeneExpBottom
  return(toExport)
}

PrepareGeneSamplesNormal <- function(GeneExpNormaldf, GeneExpTumordf, mainGene) {
  GeneExpRow <- GeneExpNormaldf[which(GeneExpTumordf[,1] == mainGene),]
  GeneTop20 <- quantile(GeneExpRow, probs = c(0.8),na.rm = TRUE)
  GeneBot20 <- quantile(GeneExpRow, probs = c(0.2),na.rm = TRUE)
  GeneExpTop <- colnames(GeneExpNormaldf)[which(GeneExpRow > GeneTop20)]
  GeneExpBottom <- colnames(GeneExpNormaldf)[which(GeneExpRow < GeneBot20)]
  toExport[[1]] <- GeneExpTop
  toExport[[2]] <- GeneExpBottom
  return(toExport)
}

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

LoadLibraries()

#dataframe containing gene data:
genedf <- data.table::fread("TCGADatasets/BRCA/probeMap_hugo_gencode_good_hg19_V24lift37_probemap")

#dataframe containing probe data for the specific tumor
probedf <- data.table::fread("probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy_Modified.csv")

#dataframe containing methylation data for samples and probes for a specific tumor
Methylationdf <- data.table::fread("TCGADatasets/BRCA/TCGA.BRCA.sampleMap_HumanMethylation45_32K.csv")
colnames(Methylationdf)[1] <- "Gene.ID"
returnListMethylation <- SeparateTumorsAndNormals(Methylationdf)
MethylationTumordf <- returnListMethylation[[1]]
MethylationNormaldf <- returnListMethylation[[2]]

GeneExpdf <- data.table::fread("TCGADatasets/BRCA/TCGA.BRCA.sampleMap_HiSeqV2.gz")
colnames(GeneExpdf)[1] <- "Gene.ID"
returnListGeneExp <- SeparateTumorsAndNormals(GeneExpdf)
GeneExpTumordf <- returnListGeneExp[[1]]
GeneExpNormaldf <- returnListGeneExp[[2]]

CNVdf <- data.table::fread("TCGADatasets/BRCA/TCGA.BRCA.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz")
colnames(CNVdf)[1] <- "Gene.ID"

Mutationdf <- data.table::fread("TCGADatasets/BRCA/mc3_gene_level_BRCA_mc3_gene_level.txt.gz")

tumorSampleOverlap <- c("Gene.ID", intersect(intersect(colnames(GeneExpTumordf)[-1], colnames(MethylationTumordf)[-1]), colnames(CNVdf)))
MethylationTumordf <- MethylationTumordf[,..tumorSampleOverlap]
GeneExpTumordf <- GeneExpTumordf[,..tumorSampleOverlap]
CNVdf <- CNVdf[,..tumorSampleOverlap]

normalSampleOverlap <- intersect(colnames(GeneExpNormaldf), colnames(MethylationNormaldf))
MethylationNormaldf <- MethylationNormaldf[,..normalSampleOverlap]
GeneExpNormaldf <- GeneExpNormaldf[,..normalSampleOverlap]

probesToStudy <- unique(intersect(unlist(Methylationdf[,1]), unlist(probedf[,1])))
genesToStudy <- unique(unlist(probedf[match(probesToStudy, unlist(probedf[,1])),2]))
if(length(which(genesToStudy == ".")) != 0) {
  genesToStudy <- genesToStudy[-(which(genesToStudy == "."))]
}

finaldfPVals <- matrix(ncol = 20)
finaldfNums <- matrix(ncol = 20)
colnames(finaldfPVals) <- c("Dataset","Gene","Location","Probe","MedianBvalTumor","MedianBvalNormal","GTMT","GTMB","GTMut","GTCG","GTCD","GBMT","GBMB","GBMut","GBCG","GBCD","GTMTnSig","GTMBnSig","GBMTnSig","GBMBnSig")
colnames(finaldfNums) <- c("Dataset","Gene","Location","Probe","MedianBvalTumor","MedianBvalNormal","GTMT","GTMB","GTMut","GTCG","GTCD","GBMT","GBMB","GBMut","GBCG","GBCD","GTMTnSig","GTMBnSig","GBMTnSig","GBMBnSig")


n <- nrow(MethylationTumordf) - 1

#FOR GENE
for(k in 1:length(genesToStudy)) {
  mainGene <- genesToStudy[k]
  geneStart <- unlist(genedf[which(genedf[,2] == mainGene),4])
  geneEnd <- unlist(genedf[which(genedf[,2] == mainGene),5])
  
  toExport <- PrepareGeneSamples(GeneExpTumordf, mainGene)
  GeneExpTop <- toExport[[1]]
  GeneExpBottom <- toExport[[2]]
  
  Mutationrow <- Mutationdf[which(Mutationdf[,1] == mainGene),-1]
  MutationSamples <- colnames(Mutationdf)[which(Mutationrow > 0) + 1]
  
  CNVrow <- CNVdf[which(CNVdf[,1] == mainGene),-1]
  CopyGainSamples <- colnames(CNVdf)[which(CNVrow > 0) + 1]
  CopyDelSamples <- colnames(CNVdf)[which(CNVrow < 0) + 1]
  
  toExport <- PrepareGeneSamplesNormal(GeneExpNormaldf, GeneExpTumordf, mainGene)
  GeneExpTopNormal <- toExport[[1]]
  GeneExpBottomNormal <- toExport[[2]]
  
  specificGeneProbes <- intersect(unlist(probedf[which(unlist(probedf[,2]) == mainGene),1]), probesToStudy)
  for(l in 1:length(specificGeneProbes)) {
    #FOR PROBE
    newrowPVals <- vector()
    newrowNums <- vector()
    
    newrowPVals <- c(newrowPVals,DATASET)
    newrowNums <- c(newrowNums,DATASET)
    
    Probe <- specificGeneProbes[l]
    rowIndex <- which(MethylationTumordf[,1] == Probe)
    
    probeStart <- unlist(probedf[which(probedf[,1] == Probe),4])
    location <- FindLocation(probeStart,geneStart,geneEnd)
    
    medianTumor <- median(as.numeric(MethylationTumordf[rowIndex,-1]))
    medianNormal <- median(as.numeric(MethylationNormaldf[rowIndex,]))

    newrowPVals <- c(newrowPVals, c(mainGene,location, Probe, medianTumor, medianNormal))
    newrowNums <- c(newrowNums, c(mainGene,location, Probe, medianTumor, medianNormal))
    
    toExport <- PrepareGeneSamples(MethylationTumordf,Probe)
    MethylationTop <- toExport[[1]]
    MethylationBottom <- toExport[[2]]
    
    x <- list(GeneExpTop,GeneExpBottom,MethylationTop,MethylationBottom,MutationSamples,CopyGainSamples,CopyDelSamples)
    for(i in 1:6) {
      for(j in (i+1):7) {
        toReturn <- Hypergeometric(x[[i]],x[[j]],ncol(MethylationTumordf))
        newrowPVals <- c(newrowPVals,toReturn[2])
        newrowNums <- c(newrowNums,toReturn[1])
      }
    }
    
    toExport <- PrepareGeneSamplesNormal(MethylationNormaldf, MethylationTumordf, Probe)
    MethylationTopNormal <- toExport[[1]]
    MethylationBottomNormal <- toExport[[2]]
    
    x <- list(GeneExpTopNormal,GeneExpBottomNormal,MethylationTopNormal,MethylationBottomNormal)
    for(i in 1:3) {
      for(j in (i+1):4) {
        toReturn <- Hypergeometric(x[[i]],x[[j]],ncol(MethylationNormaldf))
        newrowPVals <- c(newrowPVals,toReturn[2])
        newrowNums <- c(newrowNums,toReturn[1])
      }
    }
    newrowPVals <- newrowPVals[-c(7,18,27,28,33)]
    newrowNums <- newrowNums[-c(7,18,27,28,33)]
    
    newrowPVals <- newrowPVals[-c(17,18,19,20,21,22,23,24)]
    newrowNums <- newrowNums[-c(17,18,19,20,21,22,23,24)]
    
    
    finaldfPVals <- rbind(finaldfPVals,newrowPVals)
    finaldfNums <- rbind(finaldfNums,newrowNums)
  }
}

finaldfPVals <- finaldfPVals[-1,]
finaldfNums <- finaldfNums[-1,]

for(i in 7:ncol(finaldfPVals)) {
  finaldfPVals[,i] <- p.adjust(unlist(finaldfPVals[,i]), method = "fdr")
}


write.table(finaldfPVals, "GeneTableBRCAPValsMOD.csv", sep = ",", row.names = FALSE, col.names = TRUE)
write.table(finaldfNums, "GeneTableBRCANumsMOD.csv", sep = ",", row.names = FALSE, col.names = TRUE)
