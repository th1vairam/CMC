#!/bin/bash
#simulate some data
Rscript simulation.R 100 100 0.02 "data.csv"
cd ../metanetworkSynapse/
qsub -v dataFile="../CMC/data.csv",pathv="/shared/metanetworkSynapse/",sparrowZ=1,sparrow2Z=1,sparrow2ZFDR=1,lassoBIC=1,lassoAIC=1,lassoCV1se=1,lassoCVmin=1,ridgeBIC=1,ridgeAIC=1,ridgeCV1se=1,ridgeCVmin=1,genie3=1,tigress=1,numberCore=8,outputpath="/shared/testNet/" -pe orte 8 buildNet.sh
qsub -v dataFile="../CMC/data.csv",pathv="/shared/metanetworkSynapse/ARACNE/",aracne=1,aracneFull=1,correlation=1,correlationBonferroni=1,correlationFDR=1,wgcna=1,numberCore=1,outputpath="/shared/testNet/" -pe orte 8 buildNet.sh
#./pushNet.sh -a "syn4229809" -b "../CMC/code.txt" -c "../CMC/syn.txt" -defghijklmnopqv -r "none" -s "HomoZapiens" -t "bonitis" -u "blood" -x "/shared/metanetworkSynapse/pushNetworkSynapse.R"