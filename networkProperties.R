#### Function to obtain basic network properties ####
networkProperties <- function(g,  # graph for which network properties are to be calculated (in igraph format)
                              metrics = 'density'
                              # other values for metric are: c('density','diameter','avgNodeDegree','avgPathLength','avgClusteringCoefficient','centralisation')
                              ){
  
  results <- list()
  
  # Density: Shows how sparse or dense a network is, for full network the value is 1 and for empty network the value is NaN(0/0)
  if (any(metrics %in% 'density')){  
    density <- igraph::graph.density(g)
    results$density = density
  }
  
  # Diameter: Maximum distance between all the vertex in network
  if (any(metrics %in% 'diameter')){  
    diameter <- igraph::diameter(g)
    results$diameter = diameter
  }
  
  # Average node degree: Average degree of all nodes in network
  if (any(metrics %in% 'avgNodeDegree')){  
    avgNodeDegree <- mean(igraph::degree(g))
    results$avgNodeDegree = avgNodeDegree
  }
  
  # Average path length: Average distance between any two vertex in network
  if (any(metrics %in% 'avgPathLength')){  
    avgPathLength <- igraph::average.path.length(g)
    results$avgPathLength = avgPathLength
  }
  
  # Average clustering coefficeint: Gives first hand information on network clusters
  if (any(metrics %in% 'avgClusteringCoefficient')){  
    V(g)$ClusteringCoefficeint <- igraph::transitivity(g,
                                                       v=V(g),
                                                       type = 'barrat')
    avgClusteringCoefficient <- mean(V(g)$ClusteringCoefficeint,na.rm=T)
    results$avgClusteringCoefficient = avgClusteringCoefficient
  }
  
  # Centralisation: Measures how well network has a star-like topology (closer to 1 more likely the topology is star-like)
  if (any(metrics %in% 'centralisation')){  
    centralisation <- igraph::vcount(g)/(igraph::vcount(g) - 2) * (max(igraph::degree(g))/(igraph::vcount(g) - 1) - igraph::graph.density(g))
    results$centralisation = centralisation
  }
  
  return(results)
}