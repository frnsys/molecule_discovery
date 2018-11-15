Archived code of different clustering approaches. (Abandoned; no guarantee that the code still works)

Attempts:

- Training word2vec embeddings on PubMed titles and abstracts, with document representations assembled from these embeddings based on token TF-IDF scores, as described in [1], and then clustering using DBSCAN or OPTICS.
    - See `word2vec` for a modified TensorFlow word2vec skip-gram training script
    - See `cluster` for the WISDM and clustering script
    - Outcome: way too slow and memory intensive. Generating the distance matrix for DBSCAN/OPTICS would take an extremely long time. The RV coefficient, which is used for document similarity in the WISDM paper, is not a drop-in replacement for the distance metrics used in other clustering algorithms. The RV coefficient is useful because these document representations are not vectors of fixed shapes, but rather are matrices where each document may have a different shape, and standard distance metrics can't be used here.
- Generating a compound graph where nodes are compounds and an edge between compounds `A` and `B` represent co-mentions of `A` and `B` in some PubMed article or a patent. So instead of linking compounds based on the content of the articles they're mentioned in, they're linked solely on the virtue of being mentioned together in an article, under the assumption that this indicates some meaningful similarity. Once the graph is generated, a community detection algorithm (label propagation) is used to extract "clusters" from the graph.
    - See `graph.py` for this approach.
    - Right now, only using the PubMed articles rather than both articles and patents for memory reasons.
    - Originally tried using community detection algorithms (label propagation, leading eigenvector, and multilevel), but the limiting factor is that the graph is very sparse, and it seems more reasonable just to take connected components to be clusters (the number of detected communities will be at minimum the number of connected components, and generally will be more). We encode cluster labels as one-hot vectors so ideally the number of detect clusters is relatively small (on the order of 10 or 100 rather than 1000 or 10000).
    - Outcome: Still seems like it could be promising, but working with such large graphs was too slow and memory intensive.
- Training a bidirectional LSTM autoencoder on PubMed articles to learn dense representations of the articles, then cluster on those. Because these representations are vectors, standard distance metrics can be used which opens up other clustering options.
    - See `autoencoder.py` for the autoencoder script.
    - Outcome: Ended up not trying this.
