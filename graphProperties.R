## Entry point for basic graph theoretic metrics

# Load required libraries
require(igraph)
require(synapseClient)
require(data.table)
require(stringr)
require(tools)

# Needs the dev branch
library(rGithubClient)

synapseLogin()

# Source required local R files
source('./networkProperties.R')
source('./nodeRankings.R')
source('./calcNodePropertiesAndStore.R')
source('./rownameToFirstColumn.R')

# Get current commit from Github
thisFileName <- 'graphProperties.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/CMC", 
                    ref="branch", 
                    refName="master")

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('', thisFileName))

# Select input folder ids and properties to compute
Input_IDs = c('syn3526290', 'syn3526286', 'syn3526289', 'syn4549880');#
Properties = c( 'totalDegree', 'nearNeighbor', 'pageRank', 'clustCoefficient', 'betwenness', 'eigen', 'closeness');# 

# Compute metrics
Results = list()
for (property in Properties) # Evaluate loop for each property
  for (id in Input_IDs){ # Evaluate loop for each ids
    # Query to get all network files
    Networks <- synQuery(paste0('select * from file where parentId=="',id,'"'))
    
    # Create results folder
    Folder_OBJ = Folder(name = 'Node Metrics', parentId = id)
    Folder_OBJ = synStore(Folder_OBJ)
    
    # Get gene name mapping file
    GeneNames_OBJ = synGet(Networks$file.id[grep('geneName',Networks$file.name)])
    GeneNames = fread(GeneNames_OBJ@filePath, data.table=F)
    rownames(GeneNames) = GeneNames$variableNumber
    
    # Extract networks ids to evaluate
    File.IDs = Networks$file.id[!is.na(Networks$file.method) & Networks$file.method != 'sparsity']
    
    for (file.id in File.IDs){
      Results = calcNodePropertiesAndStore(file.id,
                                           property,
                                           GeneNames,
                                           parentId = Folder_OBJ$properties$id,
                                           used = GeneNames_OBJ, 
                                           executed = thisFile)
    }
  }