#### Function to obtain basic network properties ####
networkProperties <- function(g,  # graph for which network properties are to be calculated (in igraph format)
                              metrics = 'density'
                              # other values for metric are: c('density','diameter','avgNodeDegree','avgPathLength','avgClusteringCoefficient','centralisation')
                              ){
  
  results <- list()
  
  # Density: Shows how sparse or dense a network is, for full network the value is 1 and for empty network the value is NaN(0/0)
  # Centralisation: Measures how well network has a star-like topology (closer to 1 more likely the topology is star-like)
  # Average node degree: Average degree of all nodes in network
  if (any(metrics %in% c('density','centralisation', 'avgNodeDegree'))){
    results$totalDegree = igraph::graph.strength(g)
    
    nVertex = igraph::vcount(g)    
    results$maxDegree = max(results$totalDegree)
    
    results$avgNodeDegree = mean(results$totalDegree)
    results$density <- 2*sum(E(g)$weight)/(nVertex*(nVertex-1))
    results$centralisation <- nVertex/(nVertex - 2) * (results$maxDegree/(nVertex - 1) - results$density)
  }
  
  # Diameter: Maximum distance between all the vertex in network
  # Average path length: Average distance between any two vertex in network
  if (any(metrics %in% c('diameter','avgPathLength'))){  
    results$shortestPaths <- igraph::shortest.paths(g)
    results$diameter <- max(results$shortestPaths)
    results$avgPathLength <- mean(results$shortestPaths)
  }
  
  # Average clustering coefficeint: Gives first hand information on network clusters
  if (any(metrics %in% 'avgClusteringCoefficient')){  
    results$ClusteringCoefficeint <- igraph::transitivity(g,
                                                          vids=V(g),
                                                          type = 'barrat')
    results$avgClusteringCoefficient <- mean(results$ClusteringCoefficeint,na.rm=T)
  }
  
  return(results)
}