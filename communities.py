# sudo apt install -y libxml2-dev zlib1g-dev
# http://ggigraph.org/python/doc/igraph.Graph-class.html
import igraph
from scipy import sparse
from collections import defaultdict

print('Loading CID lookup...')
id2cid = [cid.strip() for cid in open('data/graph/cids.idx', 'r')]

print('Loading adjacency matrix...')
adj_mat = sparse.load_npz('data/graph/adj.npz')

# Create graph
print('Creating graph...')
vcount = max(adj_mat.shape)
print(vcount, 'nodes')
print('max edge weight', adj_mat.max())
print('mean edge weight', adj_mat.sum()/adj_mat.getnnz())
sources, targets = adj_mat.nonzero()
edgelist = list(zip(sources, targets))
print(len(edgelist), 'edges')
weights = [adj_mat[i,j] for i, j in edgelist]
print('min edge weight', min(weights))
del adj_mat
G = igraph.Graph(vcount, edgelist, edge_attrs={'weights': weights}, directed=False)

# Find communities
print('Finding communities...')
subgraphs = G.community_label_propagation(weights='weights')
labels = subgraphs.membership
comms = defaultdict(list)
for i, label in enumerate(labels):
    comms[label].append(i)

with open('data/graph/labels.dat', 'w') as f:
    f.write('\n'.join([','.join(str(id2cid[id]) for id in ids) for ids in comms.values()]))

print(len(comms), 'communities')
