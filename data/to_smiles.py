import os
from tqdm import tqdm
from glob import glob
from rdkit import Chem, RDLogger
from multiprocessing import Pool

# Silence logs
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


def to_smiles(sdf):
    path, _ = os.path.splitext(sdf)
    fname = os.path.basename(path)
    output = os.path.join('smiles', '{}.smi'.format(fname))
    if os.path.exists(output):
        return
    smiles = []
    try:
        suppl = Chem.SDMolSupplier(sdf)
        for mol in suppl:
            if mol is None: continue
            mol_id = int(mol.GetProp('PUBCHEM_COMPOUND_CID'))
            smile = Chem.MolToSmiles(mol, isomericSmiles=True)
            smiles.append('{}\t{}'.format(mol_id, smile))
        with open(output, 'w') as f:
            f.write('\n'.join(smiles))
    except RuntimeError:
        # Some SDF files couldn't be parsed properly
        pass


if __name__ == '__main__':
    if not os.path.exists('smiles'):
        os.makedirs('smiles')

    p = Pool()
    files = glob('sdf/*.sdf')
    for _ in tqdm(p.imap(to_smiles, files), total=len(files)): pass
    p.close()