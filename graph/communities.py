"""
Looks for communities in a graph (clusters).

Setup:
sudo apt install -y libxml2-dev zlib1g-dev
"""
import igraph
import numpy as np
from tqdm import tqdm
from scipy import sparse
from collections import defaultdict

LIMIT = None

print('Loading CID lookup...')
id2cid = [cid.strip() for cid in open('data/graph/cids.idx', 'r')]

print('Loading adjacency matrix...')
adj_mat = sparse.load_npz('data/graph/adj.npz')
adj_mat = adj_mat[:LIMIT,:LIMIT]
print('Compounds:', adj_mat.shape[0])
n_components, labels = sparse.csgraph.connected_components(adj_mat, directed=False, return_labels=True)
print('Connected components:', n_components)

comms = defaultdict(list)
for i, label in enumerate(labels):
    comms[label].append(i)

# Filter out isolated components
print('Filtering out isolated compounds...')
comms = {label: idxs for label, idxs in comms.items() if len(idxs) > 1}

print('After filtering:')
comm_sizes = [len(idxs) for idxs in comms.values()]
print('N compounds:', sum(comm_sizes))
print('Max:', np.max(comm_sizes))
print('Min:', np.min(comm_sizes))
print('Mean:', np.mean(comm_sizes))
print('Median:', np.median(comm_sizes))

# Get largest subgraph
# Note: the communities found are mostly limited
# by the fact that this graph is not connected;
# so we only work with the largest subgraph
# for community detection.
# Since this route was ultimately abandoned,
# future attempts could run community detection
# on each subgraph rather than just the largest one.
idxs = max(comms.values(), key=len)
x = [[idx] for idx in idxs]
adj_mat = adj_mat[x, idxs]

# ---
# Community detection
# Much slower and more memory intensive

# Create graph
print('Creating graph...')
vcount = max(adj_mat.shape)
print(vcount, 'nodes')
print('max edge weight', adj_mat.max())
print('mean edge weight', adj_mat.sum()/adj_mat.getnnz())
sources, targets = adj_mat.nonzero()
edgelist = list(zip(sources, targets))
print(len(edgelist), 'edges')
print(len(edgelist)/vcount**2, 'density')
weights = [adj_mat[i,j] for i, j in tqdm(edgelist)]
print('min edge weight', min(weights))
del adj_mat
G = igraph.Graph(vcount, edgelist, edge_attrs={'weights': weights}, directed=False)

# Find communities
print('Finding communities...')
# clusters = G.community_label_propagation(weights='weights')
clusters = G.community_leading_eigenvector(clusters=20, weights='weights')
# levels = G.community_multilevel(weights='weights', return_levels=True)
# clusters = levels[-1]

labels = clusters.membership
comms = defaultdict(list)
for i, label in enumerate(labels):
    comms[label].append(id2cid[i])

with open('data/graph/labels.dat', 'w') as f:
    f.write('\n'.join([','.join(str(cid) for cid in cids) for cids in comms.values()]))
print(len(comms), 'communities')
