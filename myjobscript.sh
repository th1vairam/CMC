#! /bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N Job1
#$ -pe orte 12
#$ -j y
#$ -l h_rt=00:30:00
mpirun -np 1 Rscript buildMpiSparrowNet.R "testData.csv" 11 "/shared/CMC/CMCSPARROW"
