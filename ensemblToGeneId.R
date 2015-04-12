ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
queryBioMart <- function(x,ensembl){
  gene1 <- NA;
  try(gene1 <- getBM(attributes='external_gene_name',filters='ensembl_gene_id',values=x,mart=ensembl),silent=T)
  if(is.null(gene1)){
    gene1<-NA
  }
  return(gene1)
}