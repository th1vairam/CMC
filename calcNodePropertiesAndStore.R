# Function to calculate node properties
calcNodePropertiesAndStore <- function(Network.id, property, GeneNames, parentId, used = NULL, executed = NULL){
  
  # load network   
  netFile <- synGet(Network.id)
  network = fread(netFile@filePath, data.table = F)
  network = network[,-(1)]
  setnames(network,'value','weight')
  network$weight = abs(network$weight)
  
  # Create a graph
  g <- graph.data.frame(network, directed=F, vertices=GeneNames[,c('variableNumber','ensemblId','geneName')])
  #   graph <- graph.adjacency(network,weighted=T,mode='undirected',diag=F)
  
  # Calculate node properties
  Metrics <- nodeRankings(g, metrics = property)
  Metrics <- rownameToFirstColumn(Metrics, 'variableNumber')
  Metrics <- merge(GeneNames, Metrics, by='variableNumber', all = T)
  
  # Write to file
  fName = stringr::str_split(netFile$properties$name,'\\.')[[1]][1]
  write.table(Metrics, 
              file = paste0(fName,'_',property,'.tsv'),
              sep = '\t', 
              col.names=T,
              row.names=F,
              quote=F)
  
  # Write metrics file to synapse
  Files <- File(paste0(fName,'_',property,'.tsv'),
                name= paste0(fName,'_',property),
                parentId = parentId)
  
  # Set Annotations
  synAnnotations <- annotations(netFile)@annotations@stringAnnotations
  synAnnotations$fileType <- 'tsv'
  synAnnotations$dataType <- 'nodeProperties'
  synAnnotations$metric <- property
  
  annotations(Files) <- synAnnotations
  
  # Store metrics in synapse
  Files <- synStore(Files,
                    used = c(used, netFile),
                    executed = executed,
                    activityName = 'Node Ranking Metrics Calculation')
  print(Files)
  
  return(Files)
}