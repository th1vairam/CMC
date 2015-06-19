makeMultiSparseNetwork <- function(sparsitySyn, networkSyn, geneSyn,uploadFolder,networkName,nNodes,fileName){
  #move to metanetworkSynapse
  library(Matrix)
  sparObj <- synGet(sparsitySyn)
  networkSyn <- synGet(networkSyn)
  geneSyn <- synGet(geneSyn)
  
  spar <- read.delim(sparObj@filePath,stringsAsFactors=F)
  network <- read.delim(networkSyn@filePath,stringsAsFactors=F)
  gene <- read.delim(geneSyn@filePath,stringsAsFactors=F)
  
  edgeListToMatrix <- function(spar,edgeList,geneName,nNodes){
    network <- matrix(0,nNodes,nNodes)
    for (i in 1:spar){
      network[edgeList[i,1],edgeList[i,2]]<-edgeList[i,3]
    }
    colnames(network) <- geneName
    rownames(network) <- geneName
    network <- Matrix(network,sparse=TRUE)
    return(network)
  }
  
  allNetworks<-sapply(spar$V2,edgeListToMatrix,geneName$ensemblId,nNodes)
  save(allNetworks,file=fileName)
  
  synObj <- File(fileName,parentId=uploadFolder)
  anno<- synGetAnnotations(networkSyn)
  
  ###UPDATE ANNOTATIONS!!
  #anno$... = ...
  
  ###Add provenance
  #act <- Activity...
  #act <- storeEntity...
  #generatedBy(synObj) <- act
  
  ###Store
  #synObj <- synStore(synObj)
  
}