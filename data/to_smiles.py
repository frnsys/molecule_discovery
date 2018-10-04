import os
from tqdm import tqdm
from glob import glob
from rdkit import Chem
from multiprocessing import Pool

# Get PubMed CIDs
cids = set()
with open('CID-PMID', 'r') as f:
    for line in f.readlines():
        line = line.strip()
        cid, pmid, _ = line.split()
        cids.add(int(cid))


def to_smiles(sdf):
    path, _ = os.path.splitext(sdf)
    fname = os.path.basename(path)

    # Skip SDF files that don't include a PubMed compound
    _, range_l, range_u = fname.split('_')
    range_l, range_u = int(range_l), int(range_u)
    if not any(cid >= range_l and cid <= range_u for cid in cids):
        return

    output = os.path.join('smiles', '{}.smi'.format(fname))
    if os.path.exists(output):
        return
    smiles = []
    suppl = Chem.SDMolSupplier(sdf)
    for mol in suppl:
        if mol is None: continue
        mol_id = int(mol.GetProp('PUBCHEM_COMPOUND_CID'))
        if mol_id not in cids:
            continue
        smile = Chem.MolToSmiles(mol, isomericSmiles=True)
        smiles.append('{}\t{}'.format(mol_id, smile))
    with open(output, 'w') as f:
        f.write('\n'.join(smiles))


if __name__ == '__main__':
    print(len(cids), 'PubMed compounds')

    if not os.path.exists('smiles'):
        os.makedirs('smiles')

    p = Pool()
    for _ in tqdm(p.imap(to_smiles, glob('sdf/*.sdf'))): pass
    p.close()