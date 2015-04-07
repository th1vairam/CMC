#!/bin/bash



#./grabCMCdata.sh
cd ../metanetworkSynapse/
qsub -v dataFile="../CMC/sczDataSVA.csv",pathv="/shared/metanetworkSynapse/",sparrow=1,aracne=0,wgcna=0,lasso=1,ridge=1,genie3=1,tigress=1,numberCore=160 -pe orte 160 buildNet.sh
qsub -v dataFile="../CMC/sczDataSVA.csv",pathv="/shared/metanetworkSynapse/ARACNE/",sparrow=0,aracne=1,wgcna=1,lasso=0,ridge=0,genie3=0,tigress=0,numberCore=160 -pe orte 160 buildNet.sh
#./pushNet.sh -p "syn3500049" -c "../CMC/code.txt" -b "../CMC/syn.txt" -sawlrgt -v "none" -o "HomoZapiens" -d "bonitis" -u "blood"
#./cleanNet.sh

