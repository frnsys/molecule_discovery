#!/bin/bash

# Download PubChem SDFs
mkdir -p sdf
cd sdf
wget -c -i ../pubchem.txt
for f in *.gz; do
    gunzip $f
done
cd ..

# Download PubMed Open Access articles
mkdir -p pubmed
cd pubmed
wget -c -i ../pubmed.txt
for f in *.tar.gz; do
    tar -xzvf $f
    rm $f
done
cd ..

# Download Compound ID -> PubMed ID associations
wget -c ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-PMID.gz
gunzip CID-PMID.gz
