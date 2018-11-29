import os
import time
import requests
import lxml.html
from collections import namedtuple
from tqdm import tqdm

Molecule = namedtuple('Molecule', ['smiles', 'is_terminal'])

def get_precursors(smiles):
    resp = requests.get('http://askcos.mit.edu/retro/target={}'.format(smiles))
    html = lxml.html.fromstring(resp.content)
    precursors = []
    for result in html.cssselect('#results tr'):
        reactants = []
        for a, term in zip(result.cssselect('.smiles a'), result.cssselect('.smiles i')):
            # href = a.attrib['href']
            smi = a.text
            terminal = 'cannot buy' not in term.text
            reactants.append((smi, terminal))
        if reactants:
            precursors.append(reactants)
    return precursors


class Node:
    def __init__(self, state, parent=None):
        self.state = state
        self.children = []
        self.parent = parent

    @property
    def is_terminal(self):
        return all(m.is_terminal for m in self.state)


def greedy_plan(target_smiles, max_depth=20):
    mol = Molecule(target_smiles, False)
    seen_states = []
    cur = Node([mol])
    path = [cur]
    for _ in tqdm(range(max_depth), desc='Search'):
        if cur.state in seen_states:
            # In a loop
            continue
        seen_states.append(cur.state)

        mols = [m for m in cur.state if not m.is_terminal]
        mol = mols[0]
        reactants = get_precursors(mol.smiles)
        new_state = list(set(cur.state) - {mol})

        if not reactants:
            # Move to end and try a diff one on next iteration
            cur.state.append(mol)
            continue

        reactants = reactants[0]
        for smi, terminal in reactants:
            new_state.append(Molecule(smi, terminal))

        cur = Node(new_state)
        path.append(cur)

        if all(s.is_terminal for s in cur.state):
            break
    else:
        print('Max depth reached')
        return path, False

    return path, True


if __name__ == '__main__':
    import json
    from multiprocessing import Pool

    def process(fname):
        with open(fname) as f:
            mols = json.load(f)
        for mol in mols:
            if 'synthesis' in mol: continue

            # Try to generate a synthesis plan
            # From base compounds -> target
            try:
                plan, complete = greedy_plan(mol['smiles'])
            except Exception as e:
                print('Exception:', e)
                time.sleep(5)
                continue
            plan = [[dict(m._asdict()) for m in n.state] for n in plan]
            mol['synthesis'] = {
                'plan': plan[::-1],
                'complete': complete
            }
        with open(fname, 'w') as f:
            json.dump(mols, f)
        time.sleep(2)

    # Get most recent batch
    out_dir = sorted(os.listdir('data/sample'))[-1]
    files = [f for f in os.listdir(out_dir) if f.endswith('.json')]
    with Pool() as p:
        for _ in tqdm(p.imap(process, files)): pass