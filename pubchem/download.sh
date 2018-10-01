#!/bin/bash

mkdir -p sdf
cd sdf
wget -i ../pubchem.txt
for f in *.gz; do
    gunzip $f
done
