#!/bin/bash

#d data file

#y file with syn ids used to run analysis
#c file with code urls used to run analysis

#s run sparrow
#a run aracne
#w run wgcna
#l run lasso
#r run ridge
#g run genie3
#t run tigress

#v whether or not sva normalization was run
#p project id of where to upload results
#b table id of where to upload network statistics

while getopts ":y:c:p:b:sawlrgtv" opt; do
  case $opt in
    y)
      echo "test y" >&2
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done