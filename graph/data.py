import json
from tqdm import tqdm
from collections import defaultdict


def stream_file(path, limit=None, skip=0):
    with open(path, 'r') as f:
        for i, line in enumerate(f):
            if i < skip:
                continue
            if limit is not None and i >= limit:
                break
            yield line


def stream_tokens(limit=None, titles_only=False):
    for line in tqdm(stream_file('text/data/pubmed.dat', limit)):
        doc = json.loads(line)
        pmid = doc['pmid']
        toks = doc['toks']
        if titles_only:
            toks = toks['title']
        else:
            toks = toks['title'] + toks['abstract']
        toks = [t.lower() for t in toks]
        yield int(pmid), toks


def load_cid2docs(articles=True, patents=False, limit=None, cids=None):
    """Group documents mentioning compounds"""
    cid2docs = defaultdict(list)
    if articles:
        for line in tqdm(stream_file('text/data/CID-PMID', limit)):
            cid, pmid, _ = line.strip().split()
            cid = int(cid)
            cid2docs[cid].append(int(pmid))
    if patents:
        for line in tqdm(stream_file('data/files/CID-Patent', limit)):
            cid, ptid = line.strip().split('\t')
            cid = int(cid)
            cid2docs[cid].append(ptid)
    if cids is not None:
        return {cid: docs for cid, docs in tqdm(cid2docs.items()) if cid in cids}
    return dict(cid2docs)


def load_cids(articles=True, patents=False, limit=None):
    cids = set()
    if articles:
        for line in tqdm(stream_file('text/data/CID-PMID', limit)):
            cid, pmid, _ = line.strip().split()
            cids.add(int(cid))
    if patents:
        for line in tqdm(stream_file('data/files/CID-Patent', limit)):
            cid, ptid = line.strip().split('\t')
            cids.add(int(cid))
    return cids
