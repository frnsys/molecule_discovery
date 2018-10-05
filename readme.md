# Setup

```
pip install -r requirements.txt
python -m spacy download en_core_web_sm
```

# Compound clustering

Attempts:

- Training word2vec embeddings on PubMed titles and abstracts, with document representations assembled from these embeddings based on token TF-IDF scores, as described in [1], and then clustering using DBSCAN or OPTICS.
    - See `data/word2vec` for a modified TensorFlow word2vec skip-gram training script
    - See `data/cluster` for the WISDM and clustering script
    - Outcome: way too slow and memory intensive. Generating the distance matrix for DBSCAN/OPTICS would take an extremely long time. The RV coefficient, which is used for document similarity in the WISDM paper, is not a drop-in replacement for the distance metrics used in other clustering algorithms. The RV coefficient is useful because these document representations are not vectors of fixed shapes, but rather are matrices where each document may have a different shape, and standard distance metrics can't be used here.
- Generating a compound graph where nodes are compounds and an edge between compounds `A` and `B` represent co-mentions of `A` and `B` in some PubMed article or a patent. So instead of linking compounds based on the content of the articles they're mentioned in, they're linked solely on the virtue of being mentioned together in an article, under the assumption that this indicates some meaningful similarity. Once the graph is generated, a community detection algorithm (label propagation) is used to extract "clusters" from the graph.
    - See `graph.py` for this approach.
    - Outcome: (still testing)
- Training a bidirectional LSTM autoencoder on PubMed articles to learn dense representations of the articles, then cluster on those. Because these representations are vectors, standard distance metrics can be used which opens up other clustering options.
    - See `autoencoder.py` for the autoencoder script.
    - Outcome: (still testing)

# Data sources

- All PubChem compounds: <ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/>
- Mapping of PubChem compound IDs to PubMed article IDs: <ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-PMID.gz>
- Mapping of PubChem compound IDs to patent IDs: <ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Patent.gz>
- PubMed Open Access articles:
    - <ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/>
    - <https://www.ncbi.nlm.nih.gov/pmc/tools/openftlist/>

# References

1. Botev, Viktor, Kaloyan Marinov, and Florian Schäfer. "Word importance-based similarity of documents metric (WISDM): Fast and scalable document similarity metric for analysis of scientific documents." Proceedings of the 6th International Workshop on Mining Scientific Publications. ACM, 2017
