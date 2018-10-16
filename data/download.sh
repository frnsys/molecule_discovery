#!/bin/bash

# Download PubChem SDFs
mkdir -p sdf
cd sdf
wget -c -i ../pubchem.txt
for f in *.gz; do
    gunzip $f
done
cd ..

# Download Compound ID -> Patent ID associations
wget -c ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Patent.gz
gunzip CID-Patent.gz
