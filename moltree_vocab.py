import os
from glob import glob
from tqdm import tqdm
from jtnn.mol_tree import MolTree, SmilesFailure
from multiprocessing import Pool
from data import load_cid2doc

# Get patent CIDs
cids = load_cid2doc(articles=False, patents=True, limit=None)
cids = set(cids.keys())
print(len(cids), 'CIDs')

def to_vocab(path):
    vocab = set()
    with open(path, 'r') as f:
        for line in f:
            cid, smiles = line.strip().split('\t')
            if int(cid) not in cids:
                continue
            try:
                mol = MolTree(smiles, skip_stereo=True)
            except SmilesFailure:
                print('SmilesFailure with CID', cid)
                continue
            for c in mol.nodes:
                vocab.add(c.smiles)
    return vocab

p = Pool()
vocab = set()
files = glob('data/smiles/*.smi')
for terms in tqdm(p.imap(to_vocab, files), total=len(files)):
    vocab = vocab.union(terms)
p.close()

print('vocab size:', len(vocab))

if not os.path.exists('data/jtnn'):
    os.makedirs('data/jtnn')

with open('data/jtnn/vocab.dat', 'w') as f:
    f.write('\n'.join(vocab))
