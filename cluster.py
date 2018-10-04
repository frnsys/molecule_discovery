import json
import math
import numpy as np
from tqdm import tqdm
from distance import dist
from itertools import product
from collections import defaultdict
from sklearn.cluster import dbscan
# from pyclustering.cluster.optics import optics;


max_n = 1000
thresh = 0.02 # 0.005

vocab = open('data/vocab.txt', 'r').read().split('\n')
embeddings = np.load('data/embeddings.npy', allow_pickle=False)
tok2id = {tok: i for i, tok in enumerate(vocab)}

idf = defaultdict(int)
tf = defaultdict(lambda: defaultdict(int))


def stream():
    with open('data/pubmed.dat', 'r') as f:
        for i, line in enumerate(f):
            if i >= max_n:
                break
            yield json.loads(line)

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
    try:
        D = np.vstack(ems)
    except ValueError:
        # Documents can't be empty
        D = np.vstack([embeddings[tok2id['UNK']]])

    pmid_idx[pmid] = len(mats)
    mats.append(D)

def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())


# Compute distance matrix
print('Computing distance matrix...')
n = len(mats)
# dist_mat = np.zeros((n, n), dtype=np.float32)
dist_mat = np.memmap('dists.dat', dtype='float32', mode='w+', shape=(n,n))

def indices(n):
    for i, j in product(range(n), range(n)):
        if j < i: yield i, j

total = ((n*n)-n)/2
for i, j in tqdm(indices(n), total=total):
    dist_mat[i, j] = dist(mats[i], mats[j])

dist_mat = symmetrize(dist_mat)
# np.save(dist_mat, 'dists.npy', allow_pickle=False)

# print('Clustering...')
_, labels = dbscan(dist_mat, eps=1., min_samples=5, metric='precomputed', n_jobs=-1)
print(len(labels))
# opt = optics(dist_mat, eps=1., minpts=5, data_type='distance_matrix')
# opt.process()
# clusters = opt.get_clusters();

# import ipdb; ipdb.set_trace()
