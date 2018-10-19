At time of writing:

- The PubChem compounds database has 134,825,000 compounds.
- 1,832,286 compounds have associated PubMed articles.
- There are 349,872,459 patent-CID associations across 6,826,637 unique patents.
- 22,626,944 compounds across the PubMed and patent CID listings.

To download PubChem compound SDF data (_warning: requires ~600GB of space_) and PubMed Open Access articles (_warning: requires ~210GB of space_):

1. `python get_urls.py` to get download URLs from PubChem.
2. `bash download.sh` to download and extract the SDF data and patent associations. Results in an `sdf` folder.

To generate derived data (_warning: this will take a several hours, days even_):

1. `python to_smiles.py` to convert SDF data for compounds into SMILES format. Results in a `smiles` folder.

To generate a sample of smiles, e.g. for quickly testing:

```bash
shuf -n 100 smiles/Compound_000000001_000025000.smi | cut -f2 > smiles_sample.txt
```