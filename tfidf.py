# msgpack-numpy==0.4.1
# msgpack-python==0.4.8

import json
import math
import numpy as np
from tqdm import tqdm
from collections import defaultdict
import msgpack_numpy as m
import msgpack

def enc_arr(arr):
    return msgpack.packb(arr, default=m.encode)

def dec_arr(enc):
    return msgpack.unpackb(enc, object_hook=m.decode)

idf = defaultdict(int)
tf = defaultdict(lambda: defaultdict(int))

thresh = 0.15

def docs():
    with open('tokens.dat', 'r') as f:
        for line in f.readlines():
            yield json.loads(line)

# TESTING
# docs = [{
#     'pmid': 0,
#     'title': 'the cat in the hat',
#     'abstract': ''
# }, {
#     'pmid': 1,
#     'title': 'the dog in the fog',
#     'abstract': ''
# }, {
#     'pmid': 2,
#     'title': 'the bat in the hat',
#     'abstract': ''
# }]


# Compute document frequencies
# and term frequencies
tok_docs = []
for doc in tqdm(docs()):
    pmid = doc['pmid']
    text = ' '.join([doc['title'], doc['abstract']])
    toks = text.split(' ')
    doc['toks'] = toks

    # Term counts
    for tok in toks:
        tf[pmid][tok] += 1

    # Document frequencies
    # and set term frequencies
    for tok in set(toks):
        idf[tok] += 1
        tf[pmid][tok] /= len(toks)
    tok_docs.append(pmid, toks)


# Compute (smoothed) inverse document frequencies
N = len(docs) + 1
for tok, count in idf.items():
    print(tok, count)
    idf[tok] = math.log(1+(N/count))


# Compute TF-IDFs
mats = {}
for pmid, toks in tok_docs:
    embeddings = []
    for tok in set(toks):
        tfidf = tf[pmid][tok] * idf[tok]
        if tfidf >= thresh:
            # TODO here lookup word2vec embeddings
            # for each of these tokens
            embeddings.append(tok)
    D = np.vstack(embeddings)

    # Not sure if this is the best way to save all these matrices
    # TODO maybe should just go ahead and compute the distance matrix here
    # without saving these representations? depends on run time
    mats[pmid] = enc_arr(D)
