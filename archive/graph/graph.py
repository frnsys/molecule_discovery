"""
Generates a compound graph
"""
import os
import numpy as np
from tqdm import tqdm
from scipy import sparse
from collections import defaultdict
from data import load_cid2docs, load_cids

m = 1.5
limit = 10000

def iter_tril(ids):
    """Iterate over the lower triangle indices pairs
    of a symmetric matrix for the given indices"""
    for i in ids:
        for j in ids:
            if j >= i: break
            yield i, j

def filter_outliers(data, m=2, key=lambda v: v):
    """Filter outlier values"""
    vals = np.array([key(v) for v in data.values()])
    mean = np.mean(vals)
    thresh = m * np.std(vals)
    return {k: v for k, v in data.items() if abs(key(v) - mean) < thresh}


if __name__ == '__main__':
    if not os.path.exists('data/graph'):
        os.makedirs('data/graph')

    # Get article-patent CIDs intersection
    cids = load_cids(articles=True, patents=False, limit=None)
    cids = cids & load_cids(articles=False, patents=True, limit=None)
    print(len(cids), 'compounds')

    # Group documents mentioning compounds
    # Not including both articles and patents b/c of memory constraints
    cid2docs = load_cid2docs(articles=False, patents=True, limit=limit, cids=cids)

    counts = [len(v) for v in cid2docs.values()]
    print(max(counts), 'most articles for an compound')
    print(np.mean(counts), 'mean articles for an compound')

    # Reduce number of compounds
    # Require at least 10 mentions
    # and filter out outliers
    # This is necessary because memory requirements became too much otherwise
    cid2docs = {cid: pids for cid, pids in cid2docs.items() if len(pids) > 10}
    cid2docs = filter_outliers(cid2docs, m=m, key=lambda v: len(v))
    print(len(cid2docs), 'compounds after filtering')

    # Save CID index for later use
    n_cids = len(cid2docs)
    id2cid = list(cid2docs.keys())
    cid2id = {cid: id for id, cid in enumerate(id2cid)}
    with open('data/graph/cids.idx', 'w') as f:
        f.write('\n'.join([str(cid) for cid in id2cid]))

    # Group compounds by their containing documents
    groups = defaultdict(list)
    for cid, pids in tqdm(cid2docs.items()):
        for pid in pids:
            id = cid2id[cid]
            groups[pid].append(id)
        groups[pid]
    counts = [len(v) for v in groups.values()]
    print(len(groups), 'articles')
    print(max(counts), 'most compounds for an article')
    print(np.mean(counts), 'mean compounds for an article')

    # Reduce number of documents
    groups = filter_outliers(groups, m=m, key=lambda v: len(v))
    print(len(groups), 'articles after filtering')

    # Create adjacency matrix
    print('Creating adjacency matrix...')
    adj_mat = sparse.lil_matrix((n_cids, n_cids), dtype=np.uint16)

    for ids in tqdm(groups.values()):
        for i, j in iter_tril(ids):
            adj_mat[i, j] += 1
    adj_mat = adj_mat + adj_mat.T - sparse.diags(adj_mat.diagonal(), dtype=np.uint16)
    print('Saving adjacency matrix...')
    sparse.save_npz('data/graph/adj.npz', adj_mat.tocsr(), compressed=True)
