"""
Model to predict ATC code for a compound.

ATC levels:

1. Anatomical main group
2. Therapeutic/pharmacological subgroup
3. Pharmacological/therapeutic/chemical subgroup
4. Pharmacological/therapeutic/chemical subgroup
5. Chemical substance

<https://www.whocc.no/atc/structure_and_principles/>
"""

import data
import random
import numpy as np

def concat(lists):
    return sum(lists, [])

def get_leaves(node):
    if isinstance(node, dict):
        return concat(map(get_leaves, node.values()))
    else:
        return list(node)

def traverse_to_depth(key, children, depth, cur_depth=1):
    if cur_depth < depth:
        return concat(traverse_to_depth(key+k, ch, depth, cur_depth=cur_depth+1)
                      for k, ch in children.items())
    else:
        return [(key, concat(map(get_leaves, children.values())))]

def split_levels(code):
    return [
        code[0],
        code[1:3],
        code[3],
        code[4],
        code[5:]
    ]

LEVEL = 3
atcs = data.load_atc_codes()
print('Coded compounds:', len(atcs))

class ATCHierarchy:
    def __init__(self, atcs):
        # Construct a hierarchy tree
        # for the ATC codes
        self._atcs = atcs
        self._hier = {}
        for cid, codes in atcs.items():
            for code in codes:
                levels = split_levels(code)
                grp = self._hier
                last = levels[-1]
                for l in levels[:-1]:
                    if l not in grp:
                        grp[l] = {}
                    grp = grp[l]
                if last not in grp:
                    grp[last] = set()
                grp[last].add(cid)

    def codes_for_level(self, level):
        return dict(concat(traverse_to_depth(k, children, level)
                           for k, children in self._hier.items()))

hier = ATCHierarchy(atcs)
codes = hier.codes_for_level(LEVEL)
print('Codes:', len(codes))
sizes = [len(c) for c in codes.values()]
print('Mean size:', np.mean(sizes))
print('Min size:', min(sizes))
print('Max size:', max(sizes))

from nnet import NNet
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from tqdm import tqdm
from sklearn import metrics

cids = set()
for vals in codes.values():
    for v in vals:
        cids.add(v)

smiles = data.load_smiles(cids)

fprints = {}
missing = set()
for cid in tqdm(cids):
    smi = smiles.get(cid)
    if smi is None:
        missing.add(cid)
        continue
    fprint = FingerprintMols.FingerprintMol(Chem.MolFromSmiles(smi),
                                            fpSize=2048, minSize=2048,
                                            minPath=1, maxPath=7,
                                            bitsPerHash=2, useHs=True,
                                            tgtDensity=0.0, branchedPaths=True,
                                            useBondOrder=True, atomInvariants=0,
                                            fromAtoms=0, atomBits=None, bitInfo=None)
    fprint = np.array([int(i) for i in fprint.ToBitString()])
    fprints[cid] = fprint

idx = {}
for i, code in enumerate(codes.keys()):
    idx[code] = i

data = []
for cid in tqdm(cids):
    x = fprints.get(cid)
    if x is None:
        continue
    c = [''.join(split_levels(code)[:LEVEL]) for code in atcs[cid]]
    y = np.zeros(len(codes))
    for code in c:
        y[idx[code]] = 1
    data.append((x, y))

random.shuffle(data)
n_testing = 200
training, testing = data[:-n_testing], data[-n_testing:]

model = NNet(n_features=2048, n_classes=len(codes), learning_rate=0.001, layers=[128,128,128])
losses = model.train(training, epochs=600)

rev = {i: code for code, i in idx.items()}

X = np.array([s[0] for s in training])
y = np.array([s[1] for s in training])
labels = model.run(X, threshold=0.3)
print('ROC AUC micro', metrics.roc_auc_score(y, labels, average='micro'))
print('ROC AUC samples', metrics.roc_auc_score(y, labels, average='samples'))
print(metrics.classification_report(y, labels))

print('---')

X = np.array([s[0] for s in testing])
y = np.array([s[1] for s in testing])
labels = model.run(X, threshold=0.3)
print('ROC AUC micro', metrics.roc_auc_score(y, labels, average='micro'))
print('ROC AUC samples', metrics.roc_auc_score(y, labels, average='samples'))
print(metrics.classification_report(y, labels))

import ipdb; ipdb.set_trace()