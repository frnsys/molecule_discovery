"""
Generate compound clusters in two steps:

1. First pass: cluster based on target receptors and action
2. Second pass: add compounds to clusters of the first pass based on structural similarity
to compounds in each cluster

Note that this is a multilabel scenario; i.e. a compound may show up in
multiple clusters if it interacts with multiple receptors.
"""
import data
import numpy as np
from tqdm import tqdm
from collections import defaultdict
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem.Fingerprints import FingerprintMols
from multiprocessing import Pool

MIN_CLUSTER_SIZE = 50

# Silence logs
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

ch = data.load_chembl()
db = data.load_drugbank()
bd = data.load_bindingdb()
stitch = data.load_stitch()
chembl2cid = data.load_chembl2cid()

smiles = {}
target_clusters = defaultdict(lambda: defaultdict(set))
for c in ch:
    id = c['chembl_id']
    tid = c['uniprot_id']
    action = c['action']
    if id is None or tid is None or action is None:
        continue
    # Try to get PubChem CID, fallback to CHEMBL ID
    id = chembl2cid.get(id, id)
    target_clusters[tid][action].add(id)
    smiles[id] = c['smiles']
for c in db:
    id = c['chembl_id']
    if id is None:
        continue
    # Try to get PubChem CID, fallback to CHEMBL ID
    id = chembl2cid.get(id, id)
    for t in c['targets']:
        for u in t['uniprot_ids']:
            tid = u['id']
            if tid is None:
                continue
            for action in t['actions']:
                if action is None:
                    continue
                target_clusters[tid][action].add(id)
    smiles[id] = c['smiles']
target_clusters = dict(target_clusters)

print('Generating fingerprints...')
fingerprints ={}
for id, smi in tqdm(smiles.items()):
    try:
        fingerprints[id] = FingerprintMols.FingerprintMol(Chem.MolFromSmiles(smi))
    except:
        continue

sizes = []
n_clusters = 0
for tid, actions in target_clusters.items():
    for action, cs in actions.items():
        sizes.append(len(cs))
        n_clusters += 1
print('Members:', sum(sizes))
print('N clusters:', n_clusters)
print('Unique compounds:', len(smiles))
print('N receptors:', len(target_clusters))
print('Mean cluster size:', np.mean(sizes))
print('Max cluster size:', max(sizes))
print('1-mem clusters:', sum(1 for s in sizes if s == 1))
print('2-mem clusters:', sum(1 for s in sizes if s == 2))

def closest_cluster(c, sim_thresh=0.5):
    try:
        mol = Chem.MolFromSmiles(c['smiles'])
        fpr = FingerprintMols.FingerprintMol(mol)
    except:
        return

    results = []
    for t in c['targets']:
        tid = t['id']

        if tid not in target_clusters:
            # results.append((c, tid, 'unknown'))
            continue

        closest = (None, None)
        for action, compounds in target_clusters[tid].items():
            sims = []
            for c_ in compounds:
                fpr_ = fingerprints.get(c_)
                sim = None
                if fpr_ is not None:
                    sim = DataStructs.FingerprintSimilarity(fpr, fpr_)
                stitch_sim = stitch[c['cid'], c_]
                if stitch_sim is not None:
                    if sim is None:
                        stitch_sim = sim
                    else:
                        sim = max(stitch_sim, sim)
                if sim is None:
                    continue
                sims.append(sim)

            # Mean similarity to cluster's compounds
            sim = sum(sims)/len(sims)
            best = closest[-1]
            if best is None or sim > best:
                closest = (action, sim)

        # We only add a compound to one cluster,
        # though we could also add them to any cluster
        # which has a mean similarity >= sim_thresh.
        if closest[-1] >= sim_thresh:
            results.append((c, tid, action))
    return results

print('Second clustering pass...')
with Pool() as p:
    for res in tqdm(p.imap(closest_cluster, bd), total=len(bd)):
        if res is None: continue
        for c, tid, action in res:
            id = c['cid']
            smiles[id] = c['smiles']
            try:
                target_clusters[tid][action].add(id)
            except KeyError:
                target_clusters[tid] = defaultdict(set)
                target_clusters[tid][action].add(id)

sizes = []
n_clusters = 0
for tid, actions in target_clusters.items():
    for action, cs in actions.items():
        sizes.append(len(cs))
        n_clusters += 1
print('Members:', sum(sizes))
print('N clusters:', n_clusters)
print('Unique compounds:', len(smiles))
print('N receptors:', len(target_clusters))
print('Mean cluster size:', np.mean(sizes))
print('Max cluster size:', max(sizes))
print('1-mem clusters:', sum(1 for s in sizes if s == 1))
print('2-mem clusters:', sum(1 for s in sizes if s == 2))

labels = defaultdict(list)
labels_idx = []
for target, actions in target_clusters.items():
    for action, compounds in actions.items():
        if len(compounds) < MIN_CLUSTER_SIZE:
            continue
        label = '{}__{}'.format(target, action)
        if label not in labels_idx:
            labels_idx.append(label)
        for id in compounds:
            labels[id].append(labels_idx.index(label))

print('Clusters after filtering:', len(labels_idx))
print('Compounds:', len(labels))

# Create outputs for training
lines = []
for id, smi in smiles.items():
    for label in labels[id]:
        lines.append('{}\t{}\t{}'.format(id, smi, label))

with open('clusters.dat', 'w') as f:
    f.write('\n'.join(lines))

with open('labels.dat', 'w') as f:
    f.write('\n'.join(labels_idx))