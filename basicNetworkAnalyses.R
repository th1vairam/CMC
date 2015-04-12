require(synapseClient)
synapseLogin()

res <- synQuery('select name,id,method from file where projectId==\'syn3455058\' and matrixType==\'sparse\'')

loadData <- function(syn){
  synMeta<-synGet(syn)
  load(synMeta@filePath)
  return(network)
}
require(metaNet)
networks <- sapply(res$file.id,loadData)
names(networks) <- res$file.method
#netStat
netStats <- sapply(networks,networkStatistics)
colnames(hubs) <- res$file.method
library(biomaRt)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

gene1 <- getBM(attributes='external_gene_name',filters='ensembl_gene_id',values=rownames(hubs)[1:5],mart=ensembl)

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

otherGene <- as.character(unlist(sapply(names(hubs[order(hubs[,'aracne'],decreasing=T),'aracne'])[1:100],queryBioMart,ensembl)))
