#!/bin/bash

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
