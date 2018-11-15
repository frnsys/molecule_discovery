# Add JTNN lib to PYTHONPATH
import os, sys
sys.path.append(os.path.join(os.getcwd(), 'jtnn'))

import json
import molvs
import torch
import rdkit
from jtnn import Vocab, JTNNVAE

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

# Load model
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
            continue
        smi = molvs.standardize_smiles(smi)
        ok.append(smi)
    samples[label] = ok

# Save generated compounds
with open('data/jtnn/compounds.json', 'w') as f:
    json.dump(samples, f)