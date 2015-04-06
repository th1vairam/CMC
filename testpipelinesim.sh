#!/bin/bash


#simulate some data

Rscript simulation.R 100 100 0.02 "data.csv"
qsub buildNet.sh -d "data.csv" -slrgt -n 16 -p "/shared/metanetworkSynapse/"
qsub buildNet.sh -d "data.csv" -aw -p "/shared/metanetworkSynapse/"
./pushNet.sh -p "syn3500049" -c code.txt -b syn.txt -sawlrgt -v "none" -o "HomoZapiens" -d "bonitis" -u "blood"