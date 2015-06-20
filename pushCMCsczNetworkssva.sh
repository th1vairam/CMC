#!/bin/sh

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N pushingcmcnetworks
#$ -e error.txt
#$ -o out.txt


#controls
#/shared/metanetworkSynapse/pushNet.sh -a "syn3526290" -b "/shared/CMC/codeControl.txt" -c "/shared/CMC/syn.txt" -defghijklmnopqv -r "None" -s "HomoSapiens" -t "Control" -u "DLPFC" -x "/shared/metanetworkSynapse/pushNetworkSynapse.R"

#cases
/shared/metanetworkSynapse/pushNet.sh -a "syn4549880" -b "/shared/CMC/codeScz.txt" -c "/shared/CMC/synSVA.txt" -defghijklmnopqv -r "SVA" -s "HomoSapiens" -t "Schizophrenia" -u "DLPFC" -x "/shared/metanetworkSynapse/pushNetworkSynapse.R"
