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
from glob import glob
from tqdm import tqdm
from collections import Counter, defaultdict

URL = 'https://www.whocc.no/atc_ddd_index/'

def get_atc_description(code):
    resp = requests.get(URL, params={'code': code})
    html = lxml.html.fromstring(resp.content)
    links = html.cssselect('#content a')[2:6]
    descs = [l.text for l in links]
    return descs

prots = data.load_protein_names()
atcs = data.load_atc_codes()

labels = [l.strip() for l in open('../data/jtnn/labels.dat', 'r')]

COMPOUNDS = defaultdict(list)
for fname in glob('../data/sample/**/*.json'):
    mols = json.load(open(fname))
    for mol in mols:
        COMPOUNDS[mol['label']].append(mol)
clusters = {}

atc_codes = set()
for label, compounds in tqdm(COMPOUNDS.items()):
    for c in compounds:
        c['image'] = c['image'].split('/')[-1]
        c['id'] = c['image'].replace('.png', '')
    label = labels[label].replace('/', '_')
    uniprot_id, action = label.split('__')
    name = prots.get(uniprot_id)
    atcs = [c['atc_code'] for c in compounds]
    clusters[label] = {
        'name': name,
        'members': compounds,
        'atc_codes': dict(Counter(atcs))
    }
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

# Copy new images/plan images
existing_images = set(os.listdir('img'))
existing_plans = set(os.listdir('plans'))
for batch in os.listdir('../data/sample'):
    batch = os.path.join('../data/sample', batch)
    images = os.path.join(batch, 'images')
    for f in os.listdir(images):
        if f not in existing_images:
            path = os.path.join(images, f)
            shutil.copyfile(path, os.path.join('img', f))
            existing_images.add(f)

    plans = os.path.join(batch, 'plans')
    for f in os.listdir(plans):
        if f not in existing_plans:
            path = os.path.join(plans, f)
            shutil.copyfile(path, os.path.join('plans', f))
            existing_plans.add(f)


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
