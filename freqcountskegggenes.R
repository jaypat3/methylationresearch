library(tidyr)
library(FedData)
library(survival)
library(ggsurvfit)

SeparateTumorsAndNormals <- function(df1) {
  x <- colnames(df1)
  y <- substr_right(x, 3)
  z <- which(y == "-11")
  
  returnList <- list()
  returnList[[1]] <- df1[,!..z]
  returnList[[2]] <- df1[,..z]
  return(returnList)
}

BLCAdf <- data.table::fread("kegggenesetanalysis/BLCA_10perc_genesetenrichment.csv")
BRCAdf <- data.table::fread("kegggenesetanalysis/BRCA_10perc_genesetenrichment.csv")
CESCdf <- data.table::fread("kegggenesetanalysis/CESC_10perc_genesetenrichment.csv")
COADdf <- data.table::fread("kegggenesetanalysis/COAD_10perc_genesetenrichment.csv")
ESCAdf <- data.table::fread("kegggenesetanalysis/ESCA_10perc_genesetenrichment.csv")
GBMdf <- data.table::fread("kegggenesetanalysis/GBM_10perc_genesetenrichment.csv")
HNSCdf <- data.table::fread("kegggenesetanalysis/HNSC_10perc_genesetenrichment.csv")
KIRCdf <- data.table::fread("kegggenesetanalysis/KIRC_10perc_genesetenrichment.csv")
KIRPdf <- data.table::fread("kegggenesetanalysis/KIRP_10perc_genesetenrichment.csv")
LGGdf <- data.table::fread("kegggenesetanalysis/LGG_10perc_genesetenrichment.csv")
LIHCdf <- data.table::fread("kegggenesetanalysis/LIHC_10perc_genesetenrichment.csv")
LUADdf <- data.table::fread("kegggenesetanalysis/LUAD_10perc_genesetenrichment.csv")
LUSCdf <- data.table::fread("kegggenesetanalysis/LUSC_10perc_genesetenrichment.csv")
OVdf <- data.table::fread("kegggenesetanalysis/OV_10perc_genesetenrichment.csv")
PAADdf <- data.table::fread("kegggenesetanalysis/PAAD_10perc_genesetenrichment.csv")
PCPGdf <- data.table::fread("kegggenesetanalysis/PCPG_10perc_genesetenrichment.csv")
PRADdf <- data.table::fread("kegggenesetanalysis/PRAD_10perc_genesetenrichment.csv")
READdf <- data.table::fread("kegggenesetanalysis/READ_10perc_genesetenrichment.csv")
SARCdf <- data.table::fread("kegggenesetanalysis/SARC_10perc_genesetenrichment.csv")
SKCMdf <- data.table::fread("kegggenesetanalysis/SKCM_10perc_genesetenrichment.csv")
STADdf <- data.table::fread("kegggenesetanalysis/STAD_10perc_genesetenrichment.csv")
TGCTdf <- data.table::fread("kegggenesetanalysis/TGCT_10perc_genesetenrichment.csv")
THCAdf <- data.table::fread("kegggenesetanalysis/THCA_10perc_genesetenrichment.csv")
THYMdf <- data.table::fread("kegggenesetanalysis/THYM_10perc_genesetenrichment.csv")
UCECdf <- data.table::fread("kegggenesetanalysis/UCEC_10perc_genesetenrichment.csv")

datasetsvec <- list(BLCAdf,BRCAdf,CESCdf,COADdf,ESCAdf,GBMdf,HNSCdf,KIRCdf,KIRPdf,LGGdf,LIHCdf,LUADdf,LUSCdf,OVdf,PAADdf,PCPGdf,PRADdf,READdf,SARCdf,SKCMdf,STADdf,TGCTdf,THCAdf,THYMdf,UCECdf)
datasetoptions <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC")

colnames(UCECdf)
columnnames <- colnames(UCECdf)[-1]
rownames <- unlist(UCECdf[,1])

newlist <- list()

for(i in 1:length(datasetsvec)) {
  newinput <- as.matrix(datasetsvec[[i]][,-1])
  colnames(newinput) <- NULL
  newinput <- ifelse(newinput <= 0.1,1,0)
  newlist[[i]] <- newinput
}

finalsum <- matrix(0,nrow(newlist[[1]]),ncol(newlist[[1]]))
for(i in 1:length(newlist)) {
  finalsum = finalsum + newlist[[i]]
}

rownames(finalsum) <- rownames
colnames(finalsum) <- columnnames

# write.table(finalsum,"frequencycountsgeneset_0.1perc.csv",sep=",",row.names=TRUE,col.names=TRUE)

#CELL_ADHESION: 32. Leukocyte: 88

indices <- vector()
for(i in 1:length(newlist)) {
  if(newlist[[i]][32,3] == 1) {
    indices <- c(indices, i)
  }
}

perc <- "10Perc.csv"
#CHANGE THIS GENESET TOO
cell_adhesion_geneset <- data.table::fread("kegggenesets/KEGG_CELL_ADHESION_MOLECULES_CAMS.v2023.2.Hs.json_geneset.csv")
cell_adhesion_geneset <- unlist(cell_adhesion_geneset)

datasetoptions <- datasetoptions[indices]

finallist <- list()
i = 1
for(whichdata in datasetoptions) {
  datasetname <- paste(paste(paste(paste(whichdata,"datasets/GeneTable",sep=""),whichdata,sep=""),"PValsMOD_",sep=""),perc,sep="")
  datasettype <- gsub(".*Table(.*)\\PVals.*", "\\1", datasetname)
  df <- data.table::fread(datasetname)
  
  gbmtset <- unlist(unique(df[df$GBMT <= 0.05 & df$GBMTnSig > 0.05 & df$GBCG > 0.05 & df$GBCD > 0.05 & df$GBMut > 0.05]$Gene))
  finallist[[i]] <- intersect(cell_adhesion_geneset,gbmtset)
  i = i + 1
}


whichdata <- vector()
pvalslog <- vector()
pvalscox <- vector()
allprobes <- vector()
allgenes <- vector()
hazard <- vector()

for(k in 1:length(datasetoptions)) {
  datasetname <- datasetoptions[k]
  
  
  genesetname <- paste(paste(paste(paste(datasetname,"datasets/GeneTable",sep=""),datasetname,sep=""),"PValsMOD_",sep=""),perc,sep="")
  genetabledf <- data.table::fread(genesetname)
  
  foldername <- paste(datasetname,"datasets/",sep="")
  setwd(foldername)
  methylationfile <- list.files()[grep("HumanMethylation",list.files(),fixed=T)]
  finalname <- paste(foldername,methylationfile,sep="")
  setwd("..")
  methylationdf <- data.table::fread(finalname)
  methylationdf <- SeparateTumorsAndNormals(methylationdf)[[1]]
  
  
  textdfname <- paste(paste("survival_",datasetname,sep=""),"_survival.txt",sep="")
  textdf <- data.table::fread(textdfname)
  
  
  
  # genetabledf <- data.table::fread("brcadatasets/GeneTableBRCAPValsMOD_10Perc.csv")
  # methylationdf <- data.table::fread("brcadatasets/TCGA.BRCA.sampleMap_HumanMethylation45_32K.csv")
  # 
  # textdf <- data.table::fread("survival_BRCA_survival.txt")
  
  #start with just one gene, one probe

  
  #finallist[[2]] for cell_adhesion, finallist[[1]] for leukocyte for brca
  for(i in 1:length(finallist[[k]])) {
    
    studygene <- finallist[[k]][i]
    subsetdf <- genetabledf[which(genetabledf$Gene == studygene),]
    relevantprobes <- unlist(subsetdf[subsetdf$GBMT < 0.05,4])
    
    for(j in 1:length(relevantprobes)) {
      probefocus <- relevantprobes[j]
      allprobes <- c(allprobes,probefocus)
      allgenes <- c(allgenes, studygene)
      whichdata <- c(whichdata,datasetoptions[k])
      # proberow <- methylationdf[methylationdf$sample == probefocus]
      
      GeneExpRow <- methylationdf[which(methylationdf[,1] == probefocus),-1]
      GeneTop20 <- quantile(GeneExpRow, probs = c(0.9),na.rm = TRUE)
      GeneBot20 <- quantile(GeneExpRow, probs = c(0.1),na.rm = TRUE)
      valuestop <- colnames(methylationdf)[which(GeneExpRow > GeneTop20)+1]
      valuesbottom <- colnames(methylationdf)[which(GeneExpRow < GeneBot20)+1]
      
      
      testdf <- textdf[unique(c(which(unlist(textdf[,1]) %in% valuestop),which(unlist(textdf[,1]) %in% valuesbottom))),]
      
      attachannotation <- vector()
      for(i in 1:nrow(testdf)) {
        if(testdf[i,1] %in% valuestop){
          attachannotation <- c(attachannotation,1)
        }
        else {
          attachannotation <- c(attachannotation,2)
        }
      }
      
      testdf <- cbind(testdf,attachannotation)
      
      
      survfit2(Surv(PFI.time,PFI) ~ attachannotation, data = testdf) %>%
        ggsurvfit() + add_confidence_interval()
      
      pvallog <- survdiff(Surv(PFI.time,PFI) ~ attachannotation, data = testdf)$p
      pvalcox <- summary(coxph(Surv(PFI.time,PFI) ~ attachannotation, data=testdf))[[7]][5]
      pvalslog <- c(pvalslog,pvallog)
      pvalscox <- c(pvalscox,pvalcox)
      hazardval <- summary(coxph(Surv(PFI.time,PFI) ~ attachannotation, data=testdf))[[7]][2]
      hazard <- c(hazard, hazardval)
    }
  }
}
pvalslogcorrected <- p.adjust(pvalslog, method = "fdr")
pvalscoxcorrected <- p.adjust(pvalscox, method = "fdr")

finaltable <- cbind(whichdata,allgenes,allprobes,pvalslog,pvalslogcorrected,pvalscox,pvalscoxcorrected,hazard)

# newdf <- data.table::fread("CellAdhesionCandidateGenes.csv")
# temprow <- colnames(newdf)
# temprow2 <- unlist(unname(newdf))
# newvec <- c(temprow, temprow2)
# 
# newlist <- list()
# for(i in 1:length(newvec)) {
#   
#   studygene <- newvec[i]
#   subsetdf <- genetabledf[which(genetabledf$Gene == studygene),]
#   # relevantprobes <- unlist(subsetdf[,4])
#   relevantprobes <- unlist(subsetdf[(subsetdf$GBMT <= 0.05 & subsetdf$GBMTnSig > 0.05 & subsetdf$GBCG > 0.05 & subsetdf$GBCD > 0.05 & subsetdf$GBMut > 0.05),4])
#   
#   for(j in 1:length(relevantprobes)) {
#     probefocus <- relevantprobes[j]
#     allprobes <- c(allprobes,probefocus)
#     allgenes <- c(allgenes, studygene)
#     proberow <- methylationdf[methylationdf$sample == probefocus]
#     
#     GeneExpRow <- methylationdf[which(methylationdf[,1] == probefocus),-1]
#     GeneTop20 <- quantile(GeneExpRow, probs = c(0.9),na.rm = TRUE)
#     GeneBot20 <- quantile(GeneExpRow, probs = c(0.1),na.rm = TRUE)
#     GeneExpTop <- colnames(methylationdf)[which(GeneExpRow > GeneTop20)+1]
#     GeneExpBottom <- colnames(methylationdf)[which(GeneExpRow < GeneBot20)+1]
#     toExport <- list()
#     toExport[[1]] <- c(studygene,probefocus)
#     toExport[[2]] <- GeneExpTop
#     toExport[[3]] <- GeneExpBottom
#     newlist <- append(newlist, toExport)
#   }
# }


