#  Entry point for basic graph theoretic metrics
require(igraph)
require(synapseClient)
synapseLogin()

# Source iGraph files
source('./networkProperties.R')
source('./nodeRankings.R')

# Query to get all network files in rda (sparse) format
Networks <- synQuery('select name,id,method from file where projectId==\'syn3455058\' and matrixType==\'sparse\'')

graphProperties <- function(Networks){
  # load network   
  netFile <- synGet(as.character(Networks['file.id']))
  load(netFile@filePath)
  
  # convert adjacency to graph
  graph <- graph.adjacency(network,weighted=T,mode='directed',diag=F)
  
  # calculate node centralities
  centralityMetrics <- nodeRankings(graph)
  centralityMetrics$geneNames <- row.names(centralityMetrics)
  centralityMetrics <- centralityMetrics[,c(dim(centralityMetrics)[2],1:(dim(centralityMetrics)[2]-1))]
  
  # write to file
  write.table(centralityMetrics,file=paste(Networks['file.method'],'Schizophrenia_SVA_nodeRankings.tsv',sep='_'),
              sep = '\t', col.names=T,row.names=F,quote=F)
  Files <- File(paste(Networks['file.method'],'Schizophrenia_SVA_nodeRankings.tsv',sep='_'),
                name= paste(Networks['file.method'],'Schizophrenia SVA nodeRankings',sep=' '),
                parentId = 'syn3455058')
  
  # set Annotations
  synAnnotations <- c(networkProperties(graph),annotations(netFile)@annotations@stringAnnotations)
  synAnnotations$fileType <- 'tsv'
  annotations(Files) <- synAnnotations
  
  # store in synapse
  File <- synStore(File,used=netFile,executed=c(),activityName='Node Ranking Metrics Calculation')
}