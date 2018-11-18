"""
Generates HTML page to explore clusters
"""

import sys; sys.path.append('../')
import data
import json
import shutil
import requests
import lxml.html
import pandas as pd
from tqdm import tqdm
from collections import Counter

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

COMPOUNDS = pd.read_csv('../data/sample/compounds.tsv', delimiter='\t')
groups = COMPOUNDS.groupby('label')
clusters = {}

atc_codes = set()
for label, compounds in tqdm(groups):
    compounds = compounds.to_dict('records')
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

shutil.copytree('../data/sample/images', 'img')

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
