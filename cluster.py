import json
import math
import numpy as np
from tqdm import tqdm
from distance import dist
from collections import defaultdict
from sklearn.cluster import dbscan


max_n = 1000
thresh = 0.005

vocab = open('data/vocab.txt', 'r').read().split('\n')
embeddings = np.load('data/embeddings.npy', allow_pickle=False)
tok2id = {tok: i for i, tok in enumerate(vocab)}

idf = defaultdict(int)
tf = defaultdict(lambda: defaultdict(int))


def stream():
    with open('data/pubmed.dat', 'r') as f:
        for i, line in enumerate(f.readlines()):
            yield json.loads(line)
            if i >= max_n:
                break

# Compute document frequencies
# and term frequencies
print('Computing DFs and TFs...')
tok_docs = []
for doc in tqdm(stream()):
    pmid = doc['pmid']
    toks = doc['toks']
    toks = [t.lower() for t in toks['title'] + toks['abstract']]

    # Term counts
    for tok in toks:
        tf[pmid][tok] += 1

    # Document frequencies
    # and set term frequencies
    for tok in set(toks):
        idf[tok] += 1
        tf[pmid][tok] /= len(toks)
    tok_docs.append((pmid, toks))


# Compute (smoothed) inverse document frequencies
print('Computing IDFs...')
N = len(tok_docs) + 1
for tok, count in tqdm(idf.items()):
    idf[tok] = math.log(1+(N/count))


# Compute TF-IDFs
# and create doc representations
print('Creating doc matrices...')
pmid_idx = {}
mats = []
for pmid, toks in tqdm(tok_docs):
    ems = []
    for tok in set(toks):
        tfidf = tf[pmid][tok] * idf[tok]
        if tfidf >= thresh:
            tid = tok2id.get(tok)
            if tid is None:
                continue
            em = embeddings[tid]
            ems.append(em)
    D = np.vstack(ems)

    # Not sure if this is the best way to save all these matrices
    # TODO maybe should just go ahead and compute the distance matrix here
    # without saving these representations? depends on run time
    # mats[pmid] = enc_arr(D)
    pmid_idx[pmid] = len(mats)
    mats.append(D)

def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())

# Compute distance matrix
print('Computing distance matrix...')
dist_mat = np.zeros((len(mats), len(mats)), dtype=np.float32)
for i, X in tqdm(enumerate(mats), total=len(mats)):
    for j, Y in enumerate(mats):
        if i == j:
            dist_mat[i, j] = 0
        else:
            d = dist(X, Y)
            dist_mat[i, j] = d
            if j > i:
                break
dist_mat = symmetrize(dist_mat)

print('Clustering...')
_, labels = dbscan(dist_mat, eps=0.1, min_samples=5, metric='precomputed', n_jobs=-1)

import ipdb; ipdb.set_trace()
