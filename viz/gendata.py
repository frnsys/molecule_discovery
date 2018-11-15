"""
Generates HTML page to explore clusters
"""

import sys; sys.path.append('../')
import os
import data
import json
import shutil
import requests
import lxml.html
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from collections import Counter, defaultdict
from rdkit.Chem import Draw
from multiprocessing import Pool

URL = 'https://www.whocc.no/atc_ddd_index/'

def get_atc_description(code):
    resp = requests.get(URL, params={'code': code})
    html = lxml.html.fromstring(resp.content)
    links = html.cssselect('#content a')[2:6]
    descs = [l.text for l in links]
    return descs

names = data.load_compound_names()
prots = data.load_protein_names()
atcs = data.load_atc_codes()

print('Loading clusters...')
clusters = defaultdict(lambda: {
    'members': []
})
labels = [l.strip() for l in open('../data/jtnn/labels.dat', 'r')]

def process_compound(line):
    id, smi, label = line.split('\t')

    atc = atcs.get(id, set())
    label = labels[int(label)].replace('/', '_')
    name = names.get(id)
    mol = Chem.MolFromSmiles(smi)
    formula = CalcMolFormula(mol)

    # Just generate all images
    # so we don't have to worry about them later
    im_fname = 'img/{}.png'.format(id)
    if not os.path.exists(im_fname):
        im = Draw.MolToImage(mol)
        im.save(im_fname)

    return label, {
        'id': id,
        'name': name,
        'formula': formula,
        'atc_codes': list([atc[:5] for atc in atc])
    }


print('Processing compounds...')
with Pool() as p:
    with open('../data/jtnn/clusters.dat', 'r') as f:
        for label, compound in tqdm(p.imap(process_compound, f)):
            clusters[label]['members'].append(compound)

print('Processing clusters...')
atc_codes = set()
for label, info in tqdm(clusters.items()):
    compounds = info['members']
    info['members'] = sorted(compounds, key=lambda c: c['name'] or 'ZZZ')

    uniprot_id, action = label.split('__')
    name = prots.get(uniprot_id)
    atcs = sum([c['atc_codes'] for c in compounds], [])
    clusters[label]['name'] = name
    clusters[label]['atc_codes'] = dict(Counter(atcs))
    for atc in atcs:
        atc_codes.add(atc)

for label, cluster in clusters.items():
    with open('clusters/{}.json'.format(label), 'w') as f:
        json.dump(cluster, f)

with open('data/clusters.json', 'w') as f:
    cluster_meta = []
    for label in sorted(list(clusters.keys())):
        cluster_meta.append({
            'label': label,
            'name': clusters[label]['name'],
            'size': len(clusters[label]['members'])
        })
    json.dump(cluster_meta, f)

# Get ATC code descriptions
try:
    with open('../data/files/atc_descs.json', 'r') as f:
        atc_descs = json.load(f)
except FileNotFoundError:
    atc_descs = {}
for code in tqdm(atc_codes):
    if code not in atc_descs:
        atc_descs[code] = get_atc_description(code)
with open('../data/files/atc_descs.json', 'w') as f:
    json.dump(atc_descs, f)
shutil.copyfile('../data/files/atc_descs.json', 'data/atc_descs.json')
