import os
import json
import molvs
from tqdm import tqdm
from glob import glob
from atc import ATCModel, code_lookup
from hashlib import md5

from rdkit import Chem, RDLogger
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from datetime import datetime

# Silence rdkit
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# Create output dir if necessary
out_dir = 'data/sample'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    os.makedirs(os.path.join(out_dir, 'images'))

# Load ATC prediction model
atc_model = ATCModel.load('data/atc')
atc_lookup = code_lookup()

# Load existing PubChem compounds to check against
pubchem = set()
for fn in tqdm(glob('data/smiles/*.smi'), desc='Loading existing compounds'):
    with open(fn, 'r') as f:
        for line in f:
            id, smi = line.strip().split('\t')
            pubchem.add(smi)


def process(fname):
    results = []
    label = int(os.path.basename(fname).replace('.dat', ''))
    with open(fname, 'r') as f:
        smis = set([l.strip() for l in f])

    ok = []
    for smi in smis:

        # Validate SMILES
        errs = molvs.validate_smiles(smi)
        if errs:
            # print('Validation error(s):', errs)
            continue

        # Standardize SMILES
        smi = molvs.standardize_smiles(smi)

        # Check if exists already
        if smi in pubchem:
            # print('Exists in PubChem')
            continue

        ok.append(smi)

    #print('Kept:', len(ok))
    atc_codes = [atc_lookup[i] for i in atc_model.predict(ok)]

    # plans = []
    # for smi in ok:
    #     plan, complete = greedy_plan(smi)
    #     plan = [[dict(m._asdict()) for m in n.state] for n in plan]
    #     plans.append((plan, complete))

    # for smi, synth_plan, atc_code in zip(ok, plans, atc_codes):
    for smi, atc_code in zip(ok, atc_codes):
        mol = Chem.MolFromSmiles(smi)
        formula = CalcMolFormula(mol)

        h = md5(smi.encode('utf8')).hexdigest()
        im = Draw.MolToImage(mol)
        im_path = os.path.join(out_dir, 'images', '{}.png'.format(h))
        im.save(im_path)

        results.append({
            'label': label,
            'smiles': smi,
            'formula': formula,
            'image': im_path,
            'atc_code': atc_code,
            'created_at': datetime.utcnow().isoformat()
        })

    # Save generated compounds
    with open(os.path.join(out_dir, '{}.json'.format(label)), 'w') as f:
        json.dump(results, f)

for _ in tqdm(map(process, glob('data/samples/*.dat')), desc='Processing samples'):
    pass