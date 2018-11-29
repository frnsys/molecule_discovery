# Add JTNN and 3N-MCTS lib to PYTHONPATH
import os, sys
sys.path.append(os.path.join(os.getcwd(), 'jtnn'))
sys.path.append(os.path.join(os.getcwd(), 'mcts'))

import json
import torch
from tqdm import tqdm
from jtnn import Vocab, JTNNVAE
from datetime import datetime

# How many compounds to generate for each class
N_SAMPLES = 20

BATCH_ID = datetime.utcnow().isoformat()

def load_jtnn():
    # Load configs & data
    model_path = 'data/jtnn/model'
    conf = json.load(open('data/jtnn/config.json', 'r'))

    vocab = [t.strip() for t in open('data/jtnn/vocab.dat')]
    vocab = Vocab(vocab)

    labels = [l.strip() for l in open('data/jtnn/labels.dat')]
    n_classes = len(labels)

    # Load JTNN VAE model
    hidden_size = conf['hidden_size']
    latent_size = conf['latent_size']
    depth = conf['depth']
    model = JTNNVAE(vocab, hidden_size, latent_size, depth, n_classes)
    saves = sorted(os.listdir(model_path))
    path = os.path.join(model_path, saves[-1])
    model.load_state_dict(torch.load(path))
    model = model.cuda()
    return model, labels


if __name__ == '__main__':
    # Create output dir if necessary
    out_dir = os.path.join('data/sample', BATCH_ID)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Load JTNN model
    model, labels = load_jtnn()

    # Generate compounds
    # torch.manual_seed(0)
    for i, label in tqdm(enumerate(labels), desc='Sampling', total=len(labels)):
        samples = []
        for _ in range(N_SAMPLES):
            smi = model.sample_prior(prob_decode=True, class_=i)
            samples.append({
                'smiles': smi
            })
        fname = os.path.join(out_dir, '{}.json'.format(i))
        with open(fname, 'w') as f:
            json.dump(samples, f)