# Add JTNN and 3N-MCTS lib to PYTHONPATH
import os, sys
sys.path.append(os.path.join(os.getcwd(), 'jtnn'))

import rdkit
import torch
from tqdm import tqdm
from torch import optim
from torch.utils.data import DataLoader
from sample import load_jtnn
from jtnn import MoleculeDataset
from atc import ATCModel, load_atc_data


lg = rdkit.RDLogger.logger()
lg.setLevel(rdkit.RDLogger.CRITICAL)


def train_atc(epochs=4000):
    """Train the ATC code prediction model"""
    atc_data, atc_idx = load_atc_data(level=3)
    if not os.path.exists('data/atc/checkpoint'):
        print('Creating new')
        model = ATCModel(n_features=2048, n_classes=len(atc_idx), learning_rate=0.001, layers=[128,128,128])
    else:
        print('Loading')
        model = ATCModel.load('data/atc')
    model.train(atc_data, epochs=epochs)
    save_path = model.save('data/atc')
    print('Model saved to', save_path)


def train_jtnn(epochs=3, batch_size=8, pretrain=True):
    """Train the JTNN CVAE model"""
    train_path = 'data/jtnn/training.dat'
    model, labels = load_jtnn()
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    scheduler = optim.lr_scheduler.ExponentialLR(optimizer, 0.9)
    scheduler.step()
    dataset = MoleculeDataset(train_path, labeled=True)
    beta = 0 if pretrain else 1.0

    for epoch in range(epochs):
        dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True,
                                num_workers=4, collate_fn=lambda x:x, drop_last=True)
        iter = tqdm(enumerate(dataloader))
        for it, batch in iter:
            for mol_tree in batch:
                mol_tree, label = mol_tree
                for node in mol_tree.nodes:
                    if node.label not in node.cands:
                        node.cands.append(node.label)
                        node.cand_mols.append(node.label_mol)

            model.zero_grad()
            loss, kl_div, wacc, tacc, sacc, dacc = model(batch, beta=beta, conditional=True)
            loss.backward()
            optimizer.step()

            iter.set_postfix(
                ep=epoch,
                kl=kl_div,
                word=wacc, # word accuracy
                topo=tacc, # topo accuracy
                assm=sacc, # assm accuracy
                steo=dacc  # steo accuracy
            )
            torch.cuda.empty_cache()

        scheduler.step()
        torch.save(model.state_dict(), 'data/jtnn/model/model.iter-' + str(epoch))


if __name__ == '__main__':
    import sys
    cmd = sys.argv[1]

    if cmd == 'jtnn_pre':
        train_jtnn(pretrain=True)
    elif cmd == 'jtnn':
        train_jtnn(pretrain=False)
    elif cmd == 'atc':
        train_atc()
    else:
        print('Unrecognized command')