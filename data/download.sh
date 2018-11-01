#!/bin/bash

# Download PubChem SDFs
python get_urls.py
mkdir -p sdf
cd sdf
wget -c -i /tmp/pubchem.txt
for f in *.gz; do
    gunzip $f
done
cd ..

mkdir -p files

# PubChem CID-ChEMBL lookup
wget -c ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-filtered.gz
gunzip CID-Synonym-filtered.gz
cat CID-Synonym-filtered | grep "\sCHEMBL" > CID-CHEMBL
rm CID-Synonym-filtered
mv CID-CHEMBL files/

# PubChem ATC codes
python get_atc.py

# BindingDB
wget -c https://www.bindingdb.org/bind/downloads/BindingDB_All_2018m9.tsv.zip
unzip BindingDB_All_2018m9.tsv.zip

# UniProt
wget -c "https://www.uniprot.org/uniprot/?query=taxonomy:%22Homo%20sapiens%20(Human)%20[9606]%22&format=tab&force=true&cols=id,entry%20name,reviewed,protein%20names,genes,organism,length&compress=yes" -O uniprot_human.tsv.gz
gunzip uniprot_human.tsv.gz
mv uniprot_human.tsv files/

# ChEMBL
wget -c ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_24_1_sqlite.tar.gz
gunzip chembl_24_1_sqlite.tar.gz

# DrugBank
echo "You will have to download the DrugBank data from <http://drugbank.ca/> and place the XML file in the `./files` folder."
