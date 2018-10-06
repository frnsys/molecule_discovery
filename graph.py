# sudo apt install -y libxml2-dev zlib1g-dev
# http://ggigraph.org/python/doc/igraph.Graph-class.html
import os
import igraph
import numpy as np
from tqdm import tqdm
from scipy import sparse
from itertools import combinations
from collections import defaultdict
from data import load_cid2doc

m = 2
limit = None
# limit = 100000

def filter_outliers(data, m=2, key=lambda v: v):
    vals = np.array([key(v) for v in data.values()])
    mean = np.mean(vals)
    thresh = m * np.std(vals)
    return {k: v for k, v in data.items() if abs(key(v) - mean) < thresh}

if not os.path.exists('data/graph'):
    os.makedirs('data/graph')

# Group documents mentioning compounds
# Not including patents b/c of memory constraints
cids = load_cid2doc(articles=True, patents=False, limit=limit)

counts = [len(v) for v in cids.values()]
print(len(cids), 'compounds')
print(max(counts), 'most articles for an compound')
print(np.mean(counts), 'mean articles for an compound')

# Reduce number of compounds
# Require at least 10 mentions
cids = {cid: pids for cid, pids in cids.items() if len(pids) > 10}
cids = filter_outliers(cids, m=m, key=lambda v: len(v))
print(len(cids), 'compounds after filtering')

# Group compounds by their containing documents
groups = defaultdict(list)
for cid, pids in tqdm(cids.items()):
    for pid in pids:
        groups[pid].append(cid)
counts = [len(v) for v in groups.values()]
print(len(groups), 'articles')
print(max(counts), 'most compounds for an article')
print(np.mean(counts), 'mean compounds for an article')

# Reduce number of documents
groups = filter_outliers(groups, m=m, key=lambda v: len(v))
print(len(groups), 'articles after filtering')

# Save CID index for later use
n_cids = len(cids)
id2cid = list(cids.keys())
cid2id = {cid: id for id, cid in enumerate(id2cid)}
with open('data/graph/cids.idx', 'w') as f:
    f.write('\n'.join([str(cid) for cid in id2cid]))

# Create adjacency matrix
print('Creating adjacency matrix...')
adj_mat = sparse.lil_matrix((n_cids, n_cids), dtype=np.uint8)
for cids in tqdm(groups.values()):
    for cid_a, cid_b in combinations(cids, 2):
        i = cid2id[cid_a]
        j = cid2id[cid_b]
        adj_mat[i, j] += 1
        adj_mat[j, i] += 1
sparse.save_npz('data/graph/adj.npz', adj_mat.tocsr(), compressed=True)

# Create graph
print('Creating graph...')
vcount = max(adj_mat.shape)
sources, targets = adj_mat.nonzero()
edgelist = list(zip(sources.tolist(), targets.tolist()))
weights = [adj_mat[i,j] for i, j in edgelist]
G = igraph.Graph(vcount, edgelist, edge_attrs={'weights': weights}, directed=False)

# Find communities
print('Finding communities...')
subgraphs = G.community_label_propagation(weights='weights')
labels = subgraphs.membership
comms = defaultdict(list)
for i, label in enumerate(labels):
    comms[label].append(i)

with open('data/graph/labels.dat', 'w') as f:
    f.write('\n'.join([','.join(str(id) for id in ids) for ids in comms.values()]))
