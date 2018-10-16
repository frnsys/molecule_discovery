At time of writing:

- There are 8,659,712 PubMed articles associated with PubChem compounds.
- There are 2,150,543 Open Access PubMed articles.
- There are 493,828 Open Access PubMed articles associated with PubChem compounds.

To download the Open Access PubMed articles (_warning: requires ~210GB of space_):

1. `python get_urls.py` to get download URLs from PubMed.
2. `bash download.sh` to download and extract the Open Access PubMed articles.

To generate derived data (_warning: this will take a several hours, days even_):

1. `python process.py` to filter PubMed articles to those with PubChem compounds and extract and tokenize their titles and abstracts. Results in a `pubmed.dat` file.
