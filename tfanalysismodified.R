library(tidyr)
library(FedData)

datasetname <- "SKCMdatasets/GeneTableSKCMPValsMOD_5Perc.csv"
pvaldf <- data.table::fread(datasetname)

datasettype <- gsub(".*Table(.*)\\PVals.*", "\\1", datasetname)
percentage <- extract_numeric(datasetname)/100
invp <- 1-percentage



SeparateTumorsAndNormals <- function(df1) {
  x <- colnames(df1)
  y <- substr_right(x, 3)
  z <- which(y == "-11")
  
  returnList <- list()
  returnList[[1]] <- df1[,!..z]
  returnList[[2]] <- df1[,..z]
  return(returnList)
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


tfupdf <- GSA::GSA.read.gmt("TFpertGEOup.gmt.txt")
TFupnames <- unlist(lapply(strsplit(tfupdf$geneset.names,"_",fixed = T),"[[",1))

tfdndf <- GSA::GSA.read.gmt("TFpertGEOdn.gmt.txt")
TFdnnames <- unlist(lapply(strsplit(tfdndf$geneset.names,"_",fixed = T),"[[",1))



genedf <- data.table::fread("SKCMdatasets/HiSeqV2")
colnames(genedf)[1] <- "Gene.ID"
returnListGeneExp <- SeparateTumorsAndNormals(genedf)
genedf <- returnListGeneExp[[1]]

methylationdf <- data.table::fread("SKCMdatasets/TCGA.SKCM.sampleMap_HumanMethylation45_31K_Cleaned.csv")
colnames(methylationdf)[1] <- "Gene.ID"
returnListMethylation <- SeparateTumorsAndNormals(methylationdf)
methylationdf <- returnListMethylation[[1]]

tumorSampleOverlap <- c("Gene.ID", intersect(colnames(genedf)[-1], colnames(methylationdf)[-1]))
methylationdf <- methylationdf[,..tumorSampleOverlap]
genedf <- genedf[,..tumorSampleOverlap]

probemap <- data.table::fread("probeMap_illuminaMethyl450_hg19_GPL16304_TCGAle_Modified.csv")

# Start by processing the TF data

human_indices <- which(grepl("HUMAN", tfupdf$geneset.names, fixed = TRUE))

type_of_change <- unlist(lapply(strsplit(tfupdf$geneset.names,"_",fixed = T),"[[",2))

KO_indices <- which("KO" == type_of_change)
KD_indices <- which("KD" == type_of_change)
OE_indices <- which("OE" == type_of_change)

combined_indices <- union(union(KO_indices,KD_indices),OE_indices)

final_tf_indices <- sort(intersect(combined_indices,human_indices))

tf_human_gene_list <- TFupnames[final_tf_indices] # equivalently dn names would work


tf_geneset_up <- list()
tf_geneset_dn <- list()
tf_change <- vector()
for (i in 1:length(final_tf_indices)) {
  index = final_tf_indices[i]
  tf_indexed_geneset_up <- c(unlist(tfupdf$genesets[index]),tfupdf$geneset.descriptions[index])
  tf_geneset_up[[i]] <- c(tf_indexed_geneset_up)
  tf_indexed_geneset_dn <- c(unlist(tfdndf$genesets[index]),tfdndf$geneset.descriptions[index])
  tf_geneset_dn[[i]] <- c(tf_indexed_geneset_dn)
  if (index %in% KO_indices) {
    tf_change <- c(tf_change, "KO")
  }
  if (index %in% KD_indices) {
    tf_change <- c(tf_change, "KD")
  }
  if (index %in% OE_indices) {
    tf_change <- c(tf_change, "OE")
  }
}

# tf_human_gene_list: TF gene
# tf_geneset_up/dn: geneset for associated TF gene

probes_in_tfgenes_indices <- which(unlist(pvaldf[,2]) %in% tf_human_gene_list)
genedf_genes <- unlist(genedf[,1])
genedf <- genedf[,-1]

# Now for the 4 category specific dissection: We start our analysis with GTMT
finalfinaldf <- data.frame(matrix(nrow=0,ncol=8))

for(s in 1:4) {
  GTMT_probes <- vector()
  GTMT_genes <- vector()
  
  dset <- ""
  if(s == 1) {
    dset <- "GTMT"
    for(index in probes_in_tfgenes_indices) {
      df_row <- pvaldf[index,]
      if (df_row$GTMT <= 0.05 & df_row$GTMTnSig > 0.05 & df_row$GTCG > 0.05 & df_row$GTCD > 0.05 & df_row$GTMut > 0.05 & df_row$GTMTnSig > 0.05) {
        GTMT_probes <- c(GTMT_probes,df_row$Probe)
        GTMT_genes <- c(GTMT_genes,df_row$Gene)
      }
    }
  }
  
  else if(s == 2) {
    dset <- "GTMB"
    for(index in probes_in_tfgenes_indices) {
      df_row <- pvaldf[index,]
      if (df_row$GTMB <= 0.05 & df_row$GTMBnSig > 0.05 & df_row$GTCG > 0.05 & df_row$GTCD > 0.05 & df_row$GTMut > 0.05 & df_row$GTMBnSig > 0.05) {
        GTMT_probes <- c(GTMT_probes,df_row$Probe)
        GTMT_genes <- c(GTMT_genes,df_row$Gene)
      }
    }
  }
  
  else if(s == 3) {
    dset <- "GBMT"
    for(index in probes_in_tfgenes_indices) {
      df_row <- pvaldf[index,]
      if (df_row$GBMT <= 0.05 & df_row$GBMTnSig > 0.05 & df_row$GBCG > 0.05 & df_row$GBCD > 0.05 & df_row$GBMut > 0.05 & df_row$GBMTnSig > 0.05) {
        GTMT_probes <- c(GTMT_probes,df_row$Probe)
        GTMT_genes <- c(GTMT_genes,df_row$Gene)
      }
    }
  }
  
  if(s == 4) {
    dset <- "GBMB"
    for(index in probes_in_tfgenes_indices) {
      df_row <- pvaldf[index,]
      if (df_row$GBMB <= 0.05 & df_row$GBMBnSig > 0.05 & df_row$GBCG > 0.05 & df_row$GBCD > 0.05 & df_row$GBMut > 0.05 & df_row$GBMBnSig > 0.05) {
        GTMT_probes <- c(GTMT_probes,df_row$Probe)
        GTMT_genes <- c(GTMT_genes,df_row$Gene)
      }
    }
  }
  
  probe_in_methylationdf_indices <- vector()
  
  for(i in 1:length(GTMT_probes)) {
    GTMT_probe <- GTMT_probes[i]
    probe_index <- which(GTMT_probes[i] == methylationdf[,1])
    if (length(probe_index) == 0) {
      probe_in_methylationdf_indices <- c(probe_in_methylationdf_indices,NA)
    }
    else {
      probe_in_methylationdf_indices <- c(probe_in_methylationdf_indices, probe_index)
    }
  }
  
  # We have now shortlisted probes which we will investigate further, ruling out the possibility of CNV and mutation.
  # From now on, our analysis is restricted to each probe. Our final table for GTMT (and other categories) will look as such:
  # Cols: Probe, Gene, status:KOKD/OE: effect on TFUp, effect on TFDn
  
  
  final_probes <- vector()
  final_genes <- vector()
  final_status <- vector()
  final_TFUppos_change <- vector()
  final_TFUpneg_change <- vector()
  final_TFDnpos_change <- vector()
  final_TFDnneg_change <- vector()
  # The end goal is to do this under a for loop; start with one probe.
  
  for(j in 1:length(GTMT_probes)) {
    probe <- GTMT_probes[j]
    gene <- GTMT_genes[j]
    
    
    index <- probe_in_methylationdf_indices[j]
    gene_row <- methylationdf[index,-1]
    top20 <- quantile(gene_row, probs = c(invp), na.rm = TRUE)
    genetop <- colnames(methylationdf)[which(gene_row > top20)+1]
    bot20 <- quantile(gene_row, probs = c(percentage), na.rm = TRUE)
    genebot <- colnames(methylationdf)[which(gene_row < bot20)+1]
    genetop_indices <- which(colnames(genedf) %in% genetop)
    genebot_indices <- which(colnames(genedf) %in% genebot)
    
    pval_gene <- vector()
    log2fc <- vector()
    for (i in 1:nrow(genedf)) {
      if(length(genetop_indices) > 0 && length(genebot_indices) > 0) {
        pval <- wilcox.test(unlist(genedf[i,..genetop_indices]),unlist(genedf[i,..genebot_indices]))$p.value
        pval_gene <- c(pval_gene, pval)
        log2fc <- c(log2fc, log2(mean(unlist(genedf[i,..genetop_indices])) / mean(unlist(genedf[i,..genebot_indices]))))
      }
      else {
        pval_gene <- c(pval_gene, 1)
        log2fc <- c(log2fc,0)
      }
    }
    
    pval_gene <- p.adjust(pval_gene, method = "fdr")
    
    poslog_genes <- genedf_genes[which(log2fc >= 0.5 & pval_gene <= 0.05)]
    neglog_genes <- genedf_genes[which(log2fc <= -0.5 & pval_gene <= 0.05)]
    
    ncoltoadd <- length(which(gene == tf_human_gene_list))
    for(i in which(gene == tf_human_gene_list)) {
      type <- tf_change[i]
      uplist <- unlist(tf_geneset_up[i])
      dnlist <- unlist(tf_geneset_dn[i])
      ##### CHANGE THESE ACCORDINGLY BASED ON MATCHUP.txt
      
      pval_hyper_uppos <- Hypergeometric(uplist,poslog_genes,nrow(genedf))[2]
      pval_hyper_upneg <- Hypergeometric(uplist,neglog_genes,nrow(genedf))[2]
      pval_hyper_dnpos <- Hypergeometric(dnlist,poslog_genes,nrow(genedf))[2]
      pval_hyper_dnneg <- Hypergeometric(dnlist,neglog_genes,nrow(genedf))[2]
      
      
      final_probes <- c(final_probes,probe)
      final_genes <- c(final_genes, gene)
      final_status <- c(final_status, type)
      final_TFUppos_change <- c(final_TFUppos_change, pval_hyper_uppos)
      final_TFUpneg_change <- c(final_TFUpneg_change, pval_hyper_upneg)
      final_TFDnpos_change <- c(final_TFDnpos_change, pval_hyper_dnpos)
      final_TFDnneg_change <- c(final_TFDnneg_change, pval_hyper_dnneg) 
    }
  }
  if(length(final_probes) > 0) {
    finaldf <- data.frame(dset, final_probes,final_genes,final_status,final_TFUppos_change,final_TFUpneg_change,final_TFDnpos_change,final_TFDnneg_change)
    # View(finaldf)
    finalfinaldf <- rbind(finalfinaldf, finaldf)
  }
}


for(i in 5:ncol(finalfinaldf)) {
  finalfinaldf[,i] <- p.adjust(unlist(finalfinaldf[,i]), method = "fdr")
}

write.table(finalfinaldf, paste(datasettype, "_TF_analysis_", percentage*100, ".csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE)
