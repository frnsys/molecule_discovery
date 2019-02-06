import os
import json
import shutil
from jinja2 import FileSystemLoader, environment

BUILD_DIR = '.build'
if os.path.isdir(BUILD_DIR):
    shutil.rmtree(BUILD_DIR)
os.mkdir(BUILD_DIR)

env = environment.Environment()
env.loader = FileSystemLoader('templates')

print('Preparing data...')
receptors = {}
clusters = json.load(open('data/clusters.json'))
atc_descs = json.load(open('data/atc_descs.json'))

receptor_descs = json.load(open('data/receptor_descs.json'))
receptor_diseases = {}
for receptor_id, data in receptor_descs.items():
    diseases = []
    for da in data:
        props = {d['label']: d.get('term') for d in da['properties']}
        disease = props['IDG Disease']
        diseases.append(disease.lower())
    receptor_diseases[receptor_id] = diseases

for clus in clusters:
    label = clus['label']
    receptor_id, action = label.split('__')
    compounds = json.load(open('data/clusters/{}.json'.format(label)))
    most_recent = sorted(compounds['members'], key=lambda c: c['created_at'], reverse=True)[0]

    if clus['name'] is None:
        short_name = receptor_id
    else:
        short_name = clus['name'].split(' (')[0]

    if receptor_id not in receptors:
        receptors[receptor_id] = {
            'name': clus['name'],
            'short_name': short_name,
            'compounds': {},
            'most_recent': most_recent,
            'diseases': receptor_diseases.get(receptor_id)
        }

    # Attach extra metadata to compounds
    for c in compounds['members']:
        c['atc_code_descs'] = atc_descs[c['atc_code']][:3]
        c['receptor'] = {
            'name': short_name,
            'id': receptor_id
        }

    # Add to receptor data
    receptors[receptor_id]['compounds'][action] = compounds

    # Update most recent compound for this receptor
    if most_recent['created_at'] > receptors[receptor_id]['most_recent']['created_at']:
        receptors[receptor_id]['most_recent'] = most_recent

print('Building pages...')

# Home page
most_recent = [r['most_recent'] for r in receptors.values()]
templ = env.get_template('index.html')
html = templ.render(compounds=most_recent[:20], receptors=receptors)
with open(os.path.join(BUILD_DIR, 'index.html'), 'w') as f:
    f.write(html)

# Receptor pages
templ = env.get_template('receptor.html')
for receptor_id, data in receptors.items():
    path = os.path.join(BUILD_DIR, receptor_id)
    os.makedirs(path, exist_ok=True)
    html = templ.render(receptor=data)
    with open(os.path.join(path, 'index.html'), 'w') as f:
        f.write(html)

# About and terms pages
for page in ['about', 'terms']:
    path = os.path.join(BUILD_DIR, page)
    os.makedirs(path, exist_ok=True)
    templ = env.get_template('{}.html'.format(page))
    html = templ.render()
    with open(os.path.join(path, 'index.html'), 'w') as f:
        f.write(html)

print('Copying files...')
# Copy other files
shutil.copytree('img', '.build/img')
shutil.copytree('plans', '.build/plans')
shutil.copyfile('style.css', '.build/style.css')
