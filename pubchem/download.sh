#!/bin/bash

mkdir -p sdf
cd sdf
wget -c -i ../pubchem.txt
for f in *.gz; do
    gunzip $f
done
cd ..
wget -c ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-PMID.gz
gunzip CID-PMID.gz
