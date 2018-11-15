# Add JTNN lib to PYTHONPATH
import os, sys
sys.path.append(os.path.join(os.getcwd(), 'jtnn'))

import json
import molvs
import torch
import rdkit
from glob import glob
from tqdm import tqdm
from jtnn import Vocab, JTNNVAE
from atc import ATCModel, code_lookup

# How many compounds to generate for each class
N_SAMPLES = 100

# Silence rdkit
lg = rdkit.RDLogger.logger()
lg.setLevel(rdkit.RDLogger.CRITICAL)

# Load configs & data
conf = json.load(open('data/jtnn/config.json', 'r'))

vocab = [t.strip() for t in open('data/jtnn/vocab.dat')]
vocab = Vocab(vocab)

labels = [l.strip() for l in open('data/jtnn/labels.dat')]
n_classes = len(labels)

# Load existing PubChem compounds to check against
pubchem = set()
for fn in tqdm(glob('data/smiles/*.smi')):
    with open(fn, 'r') as f:
        for line in f:
            pubchem.add(line.strip())

# Load ATC prediction model
atc_model = ATCModel.load('data/atc')
atc_lookup = code_lookup()

# Load JTNN VAE model
hidden_size = conf['hidden_size']
latent_size = conf['latent_size']
depth = conf['depth']
model_path = 'data/jtnn/model'
model = JTNNVAE(vocab, hidden_size, latent_size, depth, n_classes)
saves = sorted(os.listdir(model_path))
path = os.path.join(model_path, saves[-1])
model.load_state_dict(torch.load(path))
model = model.cuda()

# Generate compounds
samples = {}
torch.manual_seed(0)
for i, label in enumerate(labels):
    smis = []
    for _ in range(N_SAMPLES):
        smi = model.sample_prior(prob_decode=True, class_=i)
        smis.append(smi)
    samples[label] = smis

# Validate and standardize generated compounds
for label, smis in samples.items():
    ok = []
    for smi in smis:
        errs = molvs.validate_smiles(smi)
        if errs:
            print('Validation error(s):', errs)
            continue
        smi = molvs.standardize_smiles(smi)
        if smi in pubchem:
            print('Exists in PubChem')
            continue
        ok.append(smi)
    atc_codes = [atc_lookup[i] for i in atc_model.predict(ok)]
    samples[label] = list(zip(ok, atc_codes))

# Save generated compounds
with open('data/jtnn/compounds.json', 'w') as f:
    json.dump(samples, f)