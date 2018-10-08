At time of writing:

- the PubChem compounds database has 134,825,000 compounds.
- 1,832,286 compounds have associated PubMed articles.
- 22,626,944 compounds across the PubMed and patent CID listings.
- there are 8,659,712 PubMed articles associated with PubChem compounds.
- there are 2,150,543 Open Access PubMed articles.
- there are 493,828 Open Access PubMed articles associated with PubChem compounds.

To download PubChem compound SDF data (_warning: requires ~600GB of space_) and PubMed Open Access articles (_warning: requires ~210GB of space_):

1. `python to_urls.py` to get download URLs from PubChem and PubMed.
2. `bash download.sh` to download and extract the SDF data and PubMed articles. Results in an `sdf` folder.

To generate derived data (_warning: this will take a several hours, days even_):

1. `python to_smiles.py` to convert SDF data for compounds into SMILES format. Results in a `smiles` folder.
2. `python to_docs.py` to filter PubMed articles to those with PubChem compounds and extract and tokenize their titles and abstracts. Results in a `pubmed.dat` file.