import json
from tqdm import tqdm
from collections import defaultdict


def stream_file(path, limit=None):
    with open(path, 'r') as f:
        for i, line in enumerate(f):
            if limit is not None and i >= limit:
                break
            yield line


def stream_tokens(limit=None):
    for line in tqdm(stream_file('data/pubmed.dat', limit)):
        doc = json.loads(line)
        pmid = doc['pmid']
        toks = doc['toks']
        toks = [t.lower() for t in toks['title'] + toks['abstract']]
        yield int(pmid), toks


def load_cid2doc(articles=True, patents=False, limit=None):
    """Group documents mentioning compounds"""
    cids = defaultdict(list)
    if articles:
        for line in tqdm(stream_file('data/CID-PMID', limit)):
            cid, pmid, _ = line.strip().split()
            cids[int(cid)].append(int(pmid))
    if patents:
        for line in tqdm(stream_file('data/CID-Patent', limit)):
            cid, ptid = line.strip().split('\t')
            cids[int(cid)].append(ptid)
    return dict(cids)
