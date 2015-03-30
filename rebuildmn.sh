#!/bin/sh

cd /shared/metanetwork
git pull
cd ../
R CMD INSTALL metanetwork
cd /shared/CMC/CMCSPARROW
