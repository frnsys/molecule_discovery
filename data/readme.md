At time of writing:

- the PubChem compounds database has 134,825,000 compounds.
- 1,832,286 compounds have associated PubMed articles.
- there are 8,659,712 PubMed articles associated with compounds.
- there are 2,150,543 Open Access PubMed articles.

To download PubChem compound SDF data (_warning: requires ~600GB of space_) and PubMed Open Access articles (_warning: requires ~ of space_):

- `python get_urls.py` to get download URLs from PubChem and PubMed
- `bash download.sh` to download and extract the SDF data and PubMed articles
- `python to_smiles.py` to convert SDF data for compounds with associated PubMed articles into SMILES format (_warning: this will take a several hours_)