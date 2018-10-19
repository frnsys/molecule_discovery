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
    smiles = []
    with open(path, 'r') as f:
        for line in f:
            cid, smi = line.strip().split('\t')
            if int(cid) not in cids:
                continue
            try:
                mol = MolTree(smi, skip_stereo=True)
            except SmilesFailure:
                print('SmilesFailure with CID', cid)
                continue
            for c in mol.nodes:
                vocab.add(c.smiles)
            smiles.append('{}\t{}'.format(cid, smi))
            cids.remove(cid)
    return vocab, smiles

smiles = []
vocab = set()
files = glob('data/smiles/*.smi')
with Pool() as p:
    for terms, smis in tqdm(p.imap(to_vocab, files), total=len(files)):
        vocab = vocab.union(terms)
        smiles.extend(smis)

print('vocab size:', len(vocab))
print('smiles:', len(smiles))

if not os.path.exists('data/jtnn'):
    os.makedirs('data/jtnn')

with open('data/jtnn/vocab.dat', 'w') as f:
    f.write('\n'.join(vocab))

with open('data/jtnn/smiles.txt', 'w') as f:
    f.write('\n'.join(smiles))
