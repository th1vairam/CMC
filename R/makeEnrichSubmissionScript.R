#!usr/bin/env Rscript

# Submission Script in R
# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/home/ec2-user/Work/Github/CMC/R')

# Load libraries
library(synapseClient)
library(dplyr)

# login to synapse
synapseLogin()

# Get all files and folder
All.Files = synQuery('select * from file where projectId=="syn3455058" and fileType == "tsv" and sparsityMethod != "correlationFDR" and sparsityMethod != "wgcna"', blockSize = 100)
All.Files = All.Files$collectAll()

# Get module files
Module.Files = dplyr::filter(All.Files, is.na(file.enrichmentMethod) & file.moduleMethod == "igraph:fast_greedy" & file.sparsityMethod != "correlationFDR" & file.sparsityMethod != "wgcna")

# Get all enrichment files
Enrich.Files = dplyr::filter(All.Files, file.enrichmentMethod == "Fisher")

# Unfinished enrichment files
UEnrich.Files = Module.Files[!(paste(sapply(Module.Files$file.name, function(x){strsplit(x," ")[[1]][1]}), Module.Files$file.disease) %in%
                                 paste(sapply(Enrich.Files$file.name, function(x){strsplit(x," ")[[1]][1]}), Enrich.Files$file.disease)),]

# Make directory and write shell scripts for running these files
system('mkdir sgeEnrichSubmissions')
fp_all = file(paste('sgeEnrichSubmissions/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)
for (id in All.Files$file.id){
  fp = file (paste('/home/ec2-user/Work/Github/CMC/R/sgeEnrichSubmissions/SUB',id,sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /home/ec2-user/Work/Github/CMC/R/enrichModules.R',id), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('sgeEnrichSubmissions/allSubmissions.sh'),'a+')    
  cat(paste('qsub','-cwd','-V', paste('/home/ec2-user/Work/Github/CMC/R/sgeEnrichSubmissions/SUB',id,sep='.'),
            '-o',paste('/home/ec2-user/Work/Github/CMC/R/sgeEnrichSubmissions/SUB',id,'o',sep='.'),
            '-e',paste('/home/ec2-user/Work/Github/CMC/R/sgeEnrichSubmissions/SUB',id,'e',sep='.'),
            '-l mem=7GB'),
      file=fp_all,
      sep='\n')
  close(fp_all)
}