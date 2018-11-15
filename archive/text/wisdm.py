"""
Cluster documents based on a mixture of TF-IDF and word2vec embeddings,
as described in:

Botev, Viktor, Kaloyan Marinov, and Florian Schäfer. "Word importance-based similarity of documents metric (WISDM):
    Fast and scalable document similarity metric for analysis of scientific documents."
    Proceedings of the 6th International Workshop on Mining Scientific Publications. ACM, 2017.
"""

import numpy as np
from tqdm import tqdm
from distance import dist
from itertools import product
from sklearn.cluster import dbscan
from tfidf import tf_idf
from data import stream_tokens
# from pyclustering.cluster.optics import optics;


def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())

def indices(n):
    for i, j in product(range(n), range(n)):
        if j < i: yield i, j


if __name__ == '__main__':
    limit = 1000
    tfidf_thresh = 0.02

    vocab = open('data/vocab.txt', 'r').read().split('\n')
    embeddings = np.load('data/embeddings.npy', allow_pickle=False)
    tok2id = {tok: i for i, tok in enumerate(vocab)}

    tf, idf = tf_idf(stream_tokens(limit))

    # Compute TF-IDFs
    # and create doc representations
    print('Creating doc matrices...')
    pmid_idx = {}
    mats = []
    for pmid, toks in tqdm(stream_tokens(limit)):
        ems = []
        for tok in set(toks):
            tfidf = tf[pmid][tok] * idf[tok]
            if tfidf >= tfidf_thresh:
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

    # Compute distance matrix
    print('Computing distance matrix...')
    n = len(mats)
    # dist_mat = np.zeros((n, n), dtype=np.float32)
    dist_mat = np.memmap('dists.dat', dtype=np.float32, mode='w+', shape=(n,n))

    total = ((n*n)-n)/2
    for i, j in tqdm(indices(n), total=total):
        dist_mat[i, j] = dist(mats[i], mats[j])
    dist_mat = symmetrize(dist_mat)

    # np.save(dist_mat, 'dists.npy', allow_pickle=False)

    print('Clustering...')
    _, labels = dbscan(dist_mat, eps=1., min_samples=5, metric='precomputed', n_jobs=-1)
    print(len(labels))

    # Alternatively, OPTICS instead of DBSCAN
    # opt = optics(dist_mat, eps=1., minpts=5, data_type='distance_matrix')
    # opt.process()
    # clusters = opt.get_clusters();
