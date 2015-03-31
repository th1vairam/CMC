#!/bin/sh

Rscript makeDataMatrices.R "syn2757138" "syn2757140" "sczData.csv" "controlData.csv" "bpData.csv"
Rscript makeDataMatrices.R "syn2757147" "syn2757149" "sczDataSVA.csv" "controlDataSVA.csv" "bpDataSVA.csv"