#!/bin/bash


#simulate some data

Rscript simulation.R 100 100 0.02 "data.csv"
qsub buildNet.sh -v dataFile="../CMC/data.csv" pathv="/shared/metanetworkSynapse/" sparrow=1 aracne=0 wgcna=0 lasso=1 ridge=1 genie3=1 tigress=1 -pe orte 16
#qsub buildNet.sh -d "data.csv" -slrgt -n 16 -p "/shared/metanetworkSynapse/"
#qsub buildNet.sh -d "data.csv" -slrgt -n 16 -p "/shared/metanetworkSynapse/"
#qsub buildNet.sh -d "data.csv" -slrgt -n 16 -p "/shared/metanetworkSynapse/"
#qsub buildNet.sh -d "data.csv" -slrgt -n 16 -p "/shared/metanetworkSynapse/"
#qsub buildNet.sh -d "data.csv" -slrgt -n 16 -p "/shared/metanetworkSynapse/"
qsub ../metanetworkSynapse/buildNet.sh -v dataFile="data.csv" pathv="/shared/metanetworkSynapse/ARACNE/" sparrow=0 aracne=1 wgcna=1 lasso=0 ridge=0 genie3=0 tigress=0 numberCore=16
./pushNet.sh -p "syn3500049" -c code.txt -b syn.txt -sawlrgt -v "none" -o "HomoZapiens" -d "bonitis" -u "blood"