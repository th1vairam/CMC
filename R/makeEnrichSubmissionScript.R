#!usr/bin/env Rscript

# Submission Script in R
# Clear R console screen output
cat("\014")

# Clear R workspace
setwd('/home/ec2-user/Work/Github/CMC/R')

# Load libraries
library(synapseClient)

# login to synapse
synapseLogin()

# Get all files and folder
All.Files = synQuery('select name,id,disease,normalization from file where projectId=="syn3455058" and fileType == "tsv" and moduleMethod == "igraph:fast_greedy"')
Finished.Files = synQuery('select name,id,disease,normalization from file where projectId=="syn3455058" and fileType == "tsv" and algo == "Fisher"')

All.Files = All.Files[!(paste(tools::file_path_sans_ext(All.Files$file.name),All.Files$file.disease,All.Files$file.normalization) %in%
                          paste(sapply(Finished.Files$file.name, function(x){strsplit(x," ")[[1]][1]}), Finished.Files$file.disease, Finished.Files$file.normalization)),]

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
