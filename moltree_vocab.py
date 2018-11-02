import os
from tqdm import tqdm
from jtnn.mol_tree import MolTree, SmilesFailure
from multiprocessing import Pool


def to_vocab(line):
    vocab = set()
    cid, smi, label = line.strip().split('\t')
    try:
        mol = MolTree(smi, skip_stereo=True)
    except SmilesFailure:
        print('SmilesFailure with CID', cid)
        return
    for c in mol.nodes:
        vocab.add(c.smiles)
    return vocab

vocab = set()
with Pool() as p:
    with open('clusters.txt', 'r') as f:
        for voc in tqdm(filter(None, p.imap(to_vocab, f))):
            vocab = vocab.union(voc)
print('vocab size:', len(vocab))

if not os.path.exists('data/jtnn'):
    os.makedirs('data/jtnn')

with open('data/jtnn/vocab.dat', 'w') as f:
    f.write('\n'.join(vocab))
