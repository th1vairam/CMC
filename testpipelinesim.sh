#!/bin/bash


#simulate some data
Rscript simulation.R 100 100 0.02 "data.csv"
cd ../metanetworkSynapse/
qsub -v dataFile="../CMC/data.csv",pathv="/shared/metanetworkSynapse/",sparrow=1,aracne=0,wgcna=0,lasso=1,ridge=1,genie3=1,tigress=1,numberCore=16 -pe orte 16 buildNet.sh
qsub -v dataFile="../CMC/data.csv",pathv="/shared/metanetworkSynapse/ARACNE/",sparrow=0,aracne=1,wgcna=1,lasso=0,ridge=0,genie3=0,tigress=0,numberCore=16 -pe orte 16 buildNet.sh
./pushNet.sh -p "syn3500049" -c "../CMC/code.txt" -b "../CMC/syn.txt" -sawlrgt -v "none" -o "HomoZapiens" -d "bonitis" -u "blood"