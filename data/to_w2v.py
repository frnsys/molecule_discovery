"""Generate a word2vec training doc from pubmed.dat"""

import json
from tqdm import tqdm

with open('w2v.txt', 'w') as w:
    with open('pubmed.dat', 'r') as f:
        texts = []
        for i, line in tqdm(enumerate(f)):
            doc = json.loads(line)
            texts.append(' '.join(doc['toks']['title']))
            texts.append(' '.join(doc['toks']['abstract']))
        w.write(' '.join(texts).lower())