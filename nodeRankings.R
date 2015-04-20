#### Function to measure different node ranking metrics ####
nodeRankings <- function(g,
                         metrics = c('inDegree','outDegree','totalDegree','betwenness','clustCoefficient',
                                     'closeness','alpha','pageRank','nearNeighbor')){
  
  results <- data.frame(row.names=V(g)$name)
  
  # In degree
  if (any(metrics %in% 'inDegree')){
    results$inDegree <- degree(g,
                               v=V(g), 
                               mode = 'in', 
                               loop = T,
                               normalized = T)
  }
  
  # Out degree
  if (any(metrics %in% 'outDegree')){
    results$outDegree <- degree(g,
                                v=V(g), 
                                mode = 'out', 
                                loop = T,
                                normalized = T)
  }
  
  # Total degree
  if (any(metrics %in% 'totalDegree')){
    results$totalDegree <- degree(g,
                                  v=V(g), 
                                  mode = 'all', 
                                  loop = T,
                                  normalized = T)
  }
  
  # Betwenness centrality
  if (any(metrics %in% 'betwenness')){
    results$betweenness <- betweenness(g,
                                       v= V(g),
                                       normalized = T)
  }
  
  # Clustering coefficeint
  if (any(metrics %in% 'clustCoefficient')){
    results$clustCoefficient <- transitivity(g,
                                             v=V(g),
                                             type = 'barrat')
  }
  
  # Nearest Neighbor degree
  if (any(metrics %in% 'nearNeighbor')){
    results$nearNeighbor <- graph.knn(g,
                                      vids=V(g))$knn
  }
  
  # Page Rank
  if (any(metrics %in% 'pageRank')){
    results$pageRank <- page.rank(g,
                                  vids = V(g))$vector
  }
  
  # Alpha centrality
  if (any(metrics %in% 'alpha')){
    results$alpha <- alpha.centrality(g,
                                                nodes=V(g))
  }
  
  # Closeness centrality
  if (any(metrics %in% 'closeness')){
    results$closeness <- closeness(g,
                                   v= V(g),
                                   mode = 'all',
                                   normalized = T)
  }
}