Train word2vec embeddings from PubMed titles and abstracts.

1. Run `setup.sh` to setup your system if necessary.
2. Run `python to_w2v.py` to generate a word2vec training document from `pubmed.dat`. Results in a `w2v.txt` file.
3. Run `word2vec.py` to learn skip-gram embeddings (adapted from the [Tensorflow example code](https://github.com/tensorflow/models/blob/master/tutorials/embedding/))
