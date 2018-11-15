"""
Summarize graph communities
"""
import numpy as np
from tqdm import tqdm
from text.tfidf import tf_idf
from collections import defaultdict
from data import stream_tokens, load_cid2docs

cid2pmid = load_cid2docs(articles=True, patents=False)

with open('data/graph/cids.idx', 'r') as f:
    id2cid = f.read().split('\n')

with open('data/graph/labels.dat', 'r') as f:
    clusters = [[int(i) for i in c.split(',')] for c in f.read().split('\n')]

total = len(clusters)
sizes = [len(c) for c in clusters]
print('Clusters:', total)
print('Min size:', np.min(sizes))
print('Max size:', np.max(sizes))
print('Mean size:', np.mean(sizes))

clusters = [c for c in clusters if len(c) > 1]
print('Size=1 clusters:', total-len(clusters))
print('Size>1 clusters:', len(clusters))

sizes = [len(c) for c in clusters]
print('Mean size for size>1:', np.mean(sizes))
print('Std size for size>1:', np.std(sizes))

tf, idf = tf_idf(stream_tokens(titles_only=False))
articles = dict(stream_tokens(titles_only=False))

summaries = []
n_compounds = 0
for clus in tqdm(clusters):
    toks = defaultdict(int)
    arts = set()
    for id in clus:
        cid = int(id2cid[id])
        pmids = cid2pmid[cid]
        for pmid in pmids:
            # Articles listed in the CID-PMID lookup
            # may not be in the PubMed Open Access set,
            # in which case we won't have tokens
            for tok in articles.get(pmid, []):
                toks[tok] += tf[pmid][tok] * idf[tok]
            arts.add(pmid)
    toks = sorted(toks.items(), key=lambda x: x[1], reverse=True)
    summary = toks[:50]
    if summary:
        compounds = [int(id2cid[id]) for id in clus]
        n_compounds += len(compounds)
        summaries.append({
            'compounds': compounds,
            'summary': summary,
            'articles': ['https://www.ncbi.nlm.nih.gov/pubmed/?term={}'.format(pmid) for pmid in arts]
        })

with open('data/summaries.txt', 'w') as f:
    cluster_summaries = [
        'Final cluster count is {} across {} compounds\n\n'.format(len(summaries), n_compounds)]
    for summary in summaries:
        parts = [
            'Compounds: {}'.format(', '.join([str(cid) for cid in summary['compounds']])),
            'Terms: {}'.format(', '.join([t for t, score in summary['summary']])),
            'Articles:\n{}'.format('\n'.join(summary['articles']))
        ]
        cluster_summaries.append('\n'.join(parts))
    f.write('\n\n'.join(cluster_summaries))
