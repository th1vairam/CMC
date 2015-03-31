#! /bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N Job1
#$ -pe orte 12
mpirun -np 1 Rscript buildMpiSparrowNet.R "bpDataSVA.csv" 11 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "bpSparrowSVA.csv" "BP" "SVA" "sparrow" "syn2757147" "syn2757149" "syn2299154"
mpirun -np 1 Rscript buildMpiSparrowNet.R "controlDataSVA.csv" 11 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "controlSparrowSVA.csv" "Control" "SVA" "sparrow" "syn2757147" "syn2757149" "syn2299154"
mpirun -np 1 Rscript buildMpiSparrowNet.R "sczDataSVA.csv" 11 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "sczSparrowSVA.csv" "SCZ" "SVA" "sparrow" "syn2757147" "syn2757149" "syn2299154"
mpirun -np 1 Rscript buildMpiSparrowNet.R "bpData.csv" 11 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "bpSparrow.csv" "BP" "none" "sparrow" "syn2757138" "syn2757140" "syn2299154"
mpirun -np 1 Rscript buildMpiSparrowNet.R "controlData.csv" 11 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "controlSparrow.csv" "BP" "none" "sparrow" "syn2757138" "syn2757140" "syn2299154"
mpirun -np 1 Rscript buildMpiSparrowNet.R "sczData.csv" 11 "/shared/CMC/CMCSPARROW/"
Rscript pushNetworkSynapse.R "sczSparrow.csv" "SCZ" "none" "sparrow" "syn2757138" "syn2757140" "syn2299154"