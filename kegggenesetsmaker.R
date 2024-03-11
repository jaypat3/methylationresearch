library(rjson)

for(filename in list.files("keggfiles/")) {
  location <- paste("keggfiles/",filename,sep="")
  df <- fromJSON(file=location)
  genelist <- df[[1]]$geneSymbols
  write.table(genelist, paste(filename, "_geneset.csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE)
}
