#!/bin/bash

rm code.txt
rm syn.txt
rm synSVA.txt

#./grabCMCdata.sh
echo "syn2757147" >> synSVA.txt
echo "syn2757149" >> synSVA.txt
echo "syn2299154" >> synSVA.txt

echo "syn2299154" >> syn.txt
echo "syn2757138" >> syn.txt
echo "syn2757140" >> syn.txt

echo "https://github.com/blogsdon/CMC/blob/master/runCMCcontrols.sh" >> code.txt
echo "https://github.com/blogsdon/CMC/blob/master/grabCMCdata.sh" >> code.txt

cd ../metanetworkSynapse/
qsub -v dataFile="../CMC/controlData.csv",pathv="/shared/metanetworkSynapse/",sparrow1=1,sparrow2=1,lassoCV1se=1,lassoCV1min=1,lassoAIC=1,lassoBIC=1,ridgeCV1se=1,ridgeCV1min=1,ridgeAIC=1,ridgeBIC=1,genie3=1,tigress=1,numberCore=160 -pe orte 160 buildNet.sh
qsub -v dataFile="../CMC/controlData.csv",pathv="/shared/metanetworkSynapse/ARACNE/",aracne=1,correlation=1,numberCore=8 -pe orte 8 buildNet.sh
./pushNet.sh -a "syn3526285" -b "../CMC/code.txt" -c "../CMC/syn.txt" -defghijklmnopq -r "SVA" -s "HomoSapiens" -t "Schizophrenia" -u "DorsolateralPrefrontalCortex"
#./cleanNet.sh

