import os
from glob import glob
from tqdm import tqdm
from jtnn.mol_tree import MolTree
from multiprocessing import Pool

def to_vocab(path):
    vocab = set()
    with open(path, 'r') as f:
        for line in f:
            cid, smiles = line.strip().split('\t')
            mol = MolTree(smiles)
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
