datasetname <- "BRCA"

perc <- "10Perc.csv"

genesetname <- paste(paste(paste(paste(datasetname,"datasets/GeneTable",sep=""),datasetname,sep=""),"PValsMOD_",sep=""),perc,sep="")
df <- data.table::fread(genesetname)

genemap <- data.table::fread("brcadatasets/probeMap_hugo_gencode_good_hg19_V24lift37_probemap")
probemap <- data.table::fread("probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy_Modified.csv")

gbmtgeneset <- unlist(unique(df[df$GBMT <= 0.05 & df$GBMTnSig > 0.05 & df$GBCG > 0.05 & df$GBCD > 0.05 & df$GBMut > 0.05]$Gene))
gbmtprobeset <- unlist(unique(df[df$GBMT <= 0.05 & df$GBMTnSig > 0.05 & df$GBCG > 0.05 & df$GBCD > 0.05 & df$GBMut > 0.05]$Probe))

probemapindices <- match(gbmtprobeset, probemap$`#id`)
probestartsites <- unlist(probemap[probemapindices,4])

genestartsites <- unlist(probemap$chromStart)

collect_gene <- vector()
# for(i in 1:length(genestartsites)) {
  # collect_gene <- c(collect_gene,length(which(probestartsites %in% (genestartsites[i]-1000:genestartsites[i]))))
# }

tempvec <- vector()
for(i in 1:length(genestartsites)) {
  tempvec <- c(tempvec,(genestartsites[i]-1000):genestartsites[i])
}
