import json
import spacy
from tqdm import tqdm
from glob import glob
import pubmed_parser as pp
from multiprocessing import Pool
nlp = spacy.load('en_core_web_sm')

# Get PubMed PMIDs
pmids = set()
with open('CID-PMID', 'r') as f:
    for line in tqdm(f.readlines()):
        line = line.strip()
        cid, pmid, _ = line.split()
        pmids.add(pmid)


def extract(path):
    article = pp.parse_pubmed_xml(path)
    pmid = article['pmid']
    title = article['full_title']
    abstract = article['abstract']
    if pmid not in pmids:
        return None
    doc = {
        'pmid': pmid,
        'title': title,
        'abstract': abstract,
        'toks': {
            'title': [tok.text for tok in nlp(title)],
            'abstract': [tok.text for tok in nlp(abstract)],
        }
    }
    return doc


if __name__ == '__main__':
    p = Pool()
    save_every = 5000
    print('Globbing...')
    stream = glob('pubmed/**/*.nxml')

    with open('pubmed.dat', 'a') as f:
        docs = []
        for i, doc in enumerate(tqdm(p.imap(extract, stream))):
            if doc is None:
                continue
            docs.append(json.dumps(doc))
            if i > 0 and i % save_every == 0:
                f.write('\n'.join(docs)+'\n')
                docs = []
        if docs: f.write('\n'.join(docs))
    p.close()