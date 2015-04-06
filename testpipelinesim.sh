#!/bin/bash


#simulate some data

Rscript simulation.R 500 500 0.006 "data.csv"
qsub buildNet.sh -d "data.csv" -slrgt -n 16
qsub buildNet.sh -d "data.csv" -aw
./pushNet.sh -p "syn3500049" -c code.txt -g syn.txt -sawlrgt -v "none" -o "Homo zapiens" -d "bonitis" -u "blood"