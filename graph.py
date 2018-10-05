import numpy as np
import networkx as nx
from tqdm import tqdm
from scipy import sparse
from itertools import combinations
from collections import defaultdict
from networkx.algorithms import community

m = 2
test_max = None
# test_max = 100000


def filter_outliers(data, m=2, key=lambda v: v):
    vals = np.array([key(v) for v in data.values()])
    mean = np.mean(vals)
    thresh = m * np.std(vals)
    return {k: v for k, v in data.items() if abs(key(v) - mean) < thresh}


# Group documents mentioning compounds
cids = defaultdict(list)

with open('CID-PMID', 'r') as f:
    for i, line in enumerate(tqdm(f)):
        if test_max and i > test_max:
            break
        cid, pmid, _ = line.strip().split()
        cids[int(cid)].append('{}_{}'.format('pm', pmid))

with open('CID-Patent', 'r') as f:
    for i, line in enumerate(tqdm(f)):
        if test_max and i > test_max:
            break
        cid, ptid = line.strip().split('\t')
        cids[int(cid)].append('{}_{}'.format('pt', ptid))

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
with open('cids.idx', 'w') as f:
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
sparse.save_npz('adj.npz', adj_mat.tocsr(), compressed=True)

# Create graph
print('Creating graph...')
G = nx.from_scipy_sparse_matrix(adj_mat)

# Find communities
print('Finding communities...')
# Doesn't use weights, but much faster than the other label prop algorithm
comms = community.label_propagation_communities(G)
# comms = community.label_propagation.asyn_lpa_communities(G, weight='weight')
comms = list(comms)
with open('labels.dat', 'w') as f:
    f.write('\n'.join([','.join(str(l) for l in labels) for labels in comms]))
