"""
ATC levels:

1. Anatomical main group
2. Therapeutic/pharmacological subgroup
3. Pharmacological/therapeutic/chemical subgroup
4. Pharmacological/therapeutic/chemical subgroup
5. Chemical substance

<https://www.whocc.no/atc/structure_and_principles/>
"""

import data
from tqdm import tqdm
from .hier import ATCHierarchy
from .util import get_level
from .model import ATCModel, process_smile


def code_lookup(level=3):
    _, idx = load_atc_data(level=3)
    return {i: code for code, i in idx.items()}


def load_atc_data(level=3):
    """Load a training dataset for ATC codes,
    of (fingerprint, code ID) pairs"""
    atcs = data.load_atc_codes()
    hier = ATCHierarchy(atcs)
    codes = hier.codes_for_level(level)

    # Get CIDs for codes
    cids = set()
    for vals in codes.values():
        for v in vals:
            cids.add(v)

    # Load SMILES for CIDs
    smiles = data.load_smiles(cids)

    # Generate fingerprints for molecules
    fprints = {}
    for cid in tqdm(cids):
        smi = smiles.get(cid)
        if smi is None:
            continue

        fpr = process_smile(smi)
        fprints[cid] = fpr

    # Generate ATC code index
    # Sort to maintain consistent ordering/labels
    idx = {}
    for i, code in enumerate(sorted(codes.keys())):
        idx[code] = i

    # Assemble dataset
    training = []
    for cid in tqdm(cids):
        x = fprints.get(cid)
        if x is None:
            continue
        c = [get_level(code, level) for code in atcs[cid]]
        for code in c:
            y = idx[code]
            training.append((x, y))

    return training, idx
