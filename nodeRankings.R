#### Function to measure different node ranking metrics ####
nodeRankings <- function(g, # graph for which network properties are to be calculated (in igraph format)
                         metrics = 'totalDegree', # character vector of network properties, 
# Other options are:
    # inDegree
    # outDegree
    # totalDegree
    # betwenness
    # clustCoefficient
    # closeness
    # eigen
    # pageRank
    # nearNeighbor
                         cutoff = 2){
    
  results <- data.frame(row.names=V(g)$name)
  
  # In degree
  if (any(metrics %in% 'inDegree')){
    print('Calculating in degree...')
    results$inDegree <- igraph::graph.strength(g,
                                               vids=V(g), 
                                               mode = 'in', 
                                               loops = F)
  }
  
  # Out degree
  if (any(metrics %in% 'outDegree')){
    print('Calculating out degree...')
    results$outDegree <- igraph::strength(g,
                                          vids=V(g),
                                          mode = 'out',
                                          loops = F)
  }
  
  # Total degree: An important node is involved in a large number of interactions
  if (any(metrics %in% 'totalDegree')){
    print('Calculating total degree...')
    results$totalDegree <- igraph::strength(g,
                                            vids=V(g),
                                            mode = 'all',
                                            loops = F)
  }
  
  # Betwenness centrality: An important node will lie on a high proportion of paths between other nodes in the network.
  if (any(metrics %in% 'betwenness')){
    print('Calculating betweenness...')
    results$betweenness <- igraph::betweenness.estimate(g,
                                                        vids = V(g),
                                                        directed = F,
                                                        cutoff = cutoff)
  }
  
  # Clustering coefficeint
  if (any(metrics %in% 'clustCoefficient')){
    print('Calculating clustering coefficient...')
    results$clustCoefficient <- igraph::transitivity(g,
                                                     v=V(g),
                                                     type = 'barrat')
  }
  
  # Nearest Neighbor degree
  if (any(metrics %in% 'nearNeighbor')){
    print('Calculating k nearest neighborhood degree...')
    results$nearNeighbor <- igraph::graph.knn(g,
                                              vids=V(g))$knn
  }
  
  # Page Rank:  An important node is likely to receive more connections from other nodes.
  if (any(metrics %in% 'pageRank')){
    print('Calculating page rank...')
    results$pageRank <- igraph::page.rank(g,
                                          vids = V(g))$vector
  }
  
  # Eigen centrality: An important node is connected to important neighbours.
  if (any(metrics %in% 'eigen')){
    print('Calculating eigen centrality...')
    results$eigen <- igraph::alpha.centrality(g,
                                              nodes=V(g),
                                              sparse=T)
  }
  
  # Closeness centrality: An important node is typically “close” to, and can communicate quickly with, the other nodes in the network
  if (any(metrics %in% 'closeness')){
    print('Calculating closeness centrality ...')
    results$closeness <- igraph::closeness(g,
                                           v= V(g),
                                           mode = 'all',
                                           normalized = T)
  }
  
  return(results)
}