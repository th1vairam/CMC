#! /bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N Job1
#$ -pe orte 2
#$ -j y
#$ -l h_rt=00:30:00
mpirun -np 2 Rscript buildMpiSparrowNet.R "testData.csv" 1 "/shared/CMC/CMCSPARROW"
