To download files, run:

```
bash download.sh
```

__Warning: about 650GB of space is required.__

Once files are downloaded, you'll need to generate the derived data (i.e. convert SDF data into SMILES; __warning: this will take a several hours, days even__):

```
python to_smiles.py
```

---

# PubChem

- [Website](https://pubchem.ncbi.nlm.nih.gov/)

# UniProtKB

- Download query (for homo sapiens): <https://www.uniprot.org/uniprot/?query=taxonomy:%22Homo%20sapiens%20(Human)%20[9606]%22&format=tab&force=true&cols=id,entry%20name,reviewed,protein%20names,genes,organism,length&compress=yes>

# STITCH

- [Download](http://stitch.embl.de/download/chemical_chemical.links.v5.0.tsv.gz)
- [Readme](http://stitch.embl.de/download/README)
- [License](http://creativecommons.org/licenses/by/4.0/)

# BindingDB

- [Downloads page](https://www.bindingdb.org/bind/chemsearch/marvin/SDFdownload.jsp?all_download=yes)

# ChEMBL

- [Downloads page](http://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/)
- [Schema](http://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_24_1_schema_documentation.txt)

# DrugBank

Note: the DrugBank data is not downloaded by the `downloads.sh` script, you will need to [download it manually](https://www.drugbank.ca/releases/latest).
