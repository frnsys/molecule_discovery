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
from functools import partial
from itertools import combinations
from collections import defaultdict
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from multiprocessing import Pool

MIN_CLUSTER_SIZE = 200

# Silence logs
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


# Load data
ch = data.load_chembl()
db = data.load_drugbank()
bd = data.load_bindingdb()
stitch = data.load_stitch()
chembl2cid = data.load_chembl2cid()


def closest_cluster(c, sim_thresh=0.5):
    """Find the closest cluster to a compound `c`"""
    try:
        mol = Chem.MolFromSmiles(c['smiles'])
        fpr = FingerprintMols.FingerprintMol(mol)

        # Skip compounds with molecular weight > 500
        # according to Ghose filter for druglikeness.
        # Also helps with memory issues in JTNN
        # <https://en.wikipedia.org/wiki/Lipinski's_rule_of_five>
        wgt = Descriptors.MolWt(mol)
        if wgt > 500:
            return

        n_atoms = mol.GetNumAtoms()
        if n_atoms == 0:
            print('No-atom molecule encountered')
            raise

    except:
        return

    results = []
    for t in c['targets']:
        tid = t['id']

        # Save unknown targets
        if tid not in target_clusters:
            results.append((c, tid, 'unknown'))
            continue

        # Look for the closest cluster to this compound
        closest = (None, None)
        for action, compounds in target_clusters[tid].items():
            # Compare against each compound in the cluster
            sims = []
            for c_ in compounds:
                # Fingerprint similarity
                fpr_ = fingerprints.get(c_)
                sim = None
                if fpr_ is not None:
                    sim = DataStructs.FingerprintSimilarity(fpr, fpr_)

                # STITCH similarity
                stitch_sim = stitch[c['cid'], c_]
                if stitch_sim is not None:
                    if sim is None:
                        sim = stitch_sim
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


def summarize_clusters(clusters):
    sizes = []
    n_clusters = 0
    for tid, actions in clusters.items():
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


if __name__ == '__main__':
    smiles = {}
    target_clusters = defaultdict(lambda: defaultdict(set))

    # Get compounds from CHEMBL
    for c in ch:
        id = c['chembl_id']
        tid = c['uniprot_id']
        action = c['action']

        # Only keep compounds with an ID, a target, and an action
        if id is None or tid is None or action is None:
            continue

        # Try to get PubChem CID, fallback to CHEMBL ID
        id = chembl2cid.get(id, id)

        # Add to cluster
        target_clusters[tid][action].add(id)

        # Save SMILES for later similarity calculation
        smiles[id] = c['smiles']

    # Get compounds from BindingDB
    for c in db:
        id = c['chembl_id']
        if id is None:
            continue

        # Try to get PubChem CID, fallback to CHEMBL ID
        id = chembl2cid.get(id, id)

        # Add to cluster(s)
        for t in c['targets']:
            for u in t['uniprot_ids']:
                tid = u['id']
                if tid is None:
                    continue
                for action in t['actions']:
                    if action is None:
                        continue
                    target_clusters[tid][action].add(id)

        # Save SMILES for later similarity calculation
        smiles[id] = c['smiles']
    target_clusters = dict(target_clusters)

    # Generate fingerprints upfront, for later use
    print('Generating fingerprints...')
    fingerprints ={}
    for id, smi in tqdm(smiles.items()):
        try:
            fingerprints[id] = FingerprintMols.FingerprintMol(Chem.MolFromSmiles(smi))
        except:
            continue

    summarize_clusters(target_clusters)

    # Calculate mean similarity within
    # clusters of the same action to use as the threshold
    # for adding additional compounds
    mean_sims = []
    for tid, actions in tqdm(target_clusters.items()):
        for action, compounds in actions.items():
            # No compounds to compare
            if len(compounds) < 2:
                continue

            # Calculate intracluster similarites
            sims = []
            for c, c_ in combinations(compounds, 2):
                fpr = fingerprints.get(c)
                fpr_ = fingerprints.get(c_)
                if fpr is not None and fpr_ is not None:
                    sim = DataStructs.FingerprintSimilarity(fpr, fpr_)
                    sims.append(sim)
            mean_sims.append(sum(sims)/len(sims))

    # Use 25th percentile as threshold
    threshold = np.percentile(mean_sims, 25)
    clos_clus = partial(closest_cluster, sim_thresh=threshold)

    # Now run second clustering pass, clustering in
    # compounds which did not have targets and actions
    # explicitly labeled
    print('Second clustering pass...')
    with Pool() as p:
        updates = []
        for res in tqdm(p.imap(clos_clus, bd), total=len(bd)):
            if res is None: continue
            updates += res

    # Add the compounds to their computed clusters
    for c, tid, action in updates:
        id = c['cid']
        smiles[id] = c['smiles']
        try:
            target_clusters[tid][action].add(id)
        except KeyError:
            target_clusters[tid] = defaultdict(set)
            target_clusters[tid][action].add(id)

    summarize_clusters(target_clusters)

    # Filter to clusters of a minimum size
    # (to reduce number of clusters to a more manageable
    # number and to cut down on noise)
    labels_idx = []
    labels = defaultdict(list)
    for target, actions in target_clusters.items():
        for action, compounds in actions.items():
            if len(compounds) < MIN_CLUSTER_SIZE:
                continue

            # Come up with informative labels for each cluster
            label = '{}__{}'.format(target, action)
            if label not in labels_idx:
                labels_idx.append(label)

            # Add compounds to their labeled clusters
            for id in compounds:
                labels[id].append(labels_idx.index(label))

    print('Clusters after filtering:', len(labels_idx))
    print('Compounds:', len(labels))

    # Create outputs for training
    lines = []
    for id, smi in smiles.items():
        for label in labels[id]:
            lines.append('{}\t{}\t{}'.format(id, smi, label))

    with open('data/jtnn/clusters.dat', 'w') as f:
        f.write('\n'.join(lines))

    with open('data/jtnn/labels.dat', 'w') as f:
        f.write('\n'.join(labels_idx))