#! /bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N Job1
#$ -pe orte 288
mpirun -np 1 Rscript buildMpiSparrowNet.R "bpDataSVA.csv" 287 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "bpSparrowSVA.csv" "BP" "SVA" "sparrow" "syn2757147" "syn2757149" "syn2299154"
rm bpSparrowSVA.csv
rm result.rda
mpirun -np 1 Rscript buildMpiSparrowNet.R "controlDataSVA.csv" 287 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "controlSparrowSVA.csv" "Control" "SVA" "sparrow" "syn2757147" "syn2757149" "syn2299154"
rm controlSparrowSVA.csv
rm result.rda
mpirun -np 1 Rscript buildMpiSparrowNet.R "sczDataSVA.csv" 287 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "sczSparrowSVA.csv" "SCZ" "SVA" "sparrow" "syn2757147" "syn2757149" "syn2299154"
rm sczSparrowSVA.csv
rm result.rda
mpirun -np 1 Rscript buildMpiSparrowNet.R "bpData.csv" 287 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "bpSparrow.csv" "BP" "none" "sparrow" "syn2757138" "syn2757140" "syn2299154"
rm bpSparrow.csv
rm result.rda
mpirun -np 1 Rscript buildMpiSparrowNet.R "controlData.csv" 287 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "controlSparrow.csv" "BP" "none" "sparrow" "syn2757138" "syn2757140" "syn2299154"
rm controlSparrow.csv
rm result.rda
mpirun -np 1 Rscript buildMpiSparrowNet.R "sczData.csv" 287 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "sczSparrow.csv" "SCZ" "none" "sparrow" "syn2757138" "syn2757140" "syn2299154"
rm sczSparrow.csv
rm result.rda