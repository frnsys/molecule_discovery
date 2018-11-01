"""
Generates HTML page to explore clusters
"""

import data
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from collections import defaultdict
from jinja2 import Environment, BaseLoader

names = data.load_compound_names()
prots = data.load_protein_names()
atcs = data.load_atc_codes()

clusters = defaultdict(list)
with open('clusters.txt', 'r') as f:
    for line in tqdm(f):
        id, label, smi = line.split('\t')
        name = names.get(id)
        atc = atcs.get(id)
        if name is None:
            name = CalcMolFormula(Chem.MolFromSmiles(smi))
        clusters[label].append({
            'name': name,
            'atc_code': atc
        })

cluster_meta = {}
for label in clusters.keys():
    uniprot_id, action = label.split('__')
    name = prots.get(uniprot_id)
    cluster_meta[label] = {
        'name': name
    }

tmpl = '''
<html>
<head>
<style>
html, body {
    font-family: 'Arial', sans-serif;
}
ul {
    list-style-type: none;
}
li {
    font-size: 0.75em;
    padding: 0.2em;
    margin: 0 0 0.1em 0;
    display: inline-block;
    background: #ff6c6c;
    border-radius: 0.2em;
    border: 1px solid #b34b4b;
    cursor: help;
}
.meta {
    position: absolute;
    background: #e6e237;
    color: #000;
    padding: 0.5em;
    border: 1px solid #141316;
    display: none;
}
li:hover .meta {
    display: block;
}
</style>
</head>
<body>
    {% for label, mems in clusters.items() %}
        <h2>{{ label }}</h2>
        <h3>{{ meta[label]['name'] }}</h3>
        <ul>
            {% for mem in mems %}
                <li>
                    {{ mem.name }}
                    <div class="meta">ATC: {{ mem.atc_code }}</div>
                </li>
            {% endfor %}
        </ul>
    {% endfor %}
</body>
</html>
'''

tmpl = Environment(loader=BaseLoader()).from_string(source=tmpl)
html = tmpl.render(clusters=clusters, meta=cluster_meta)

with open('clusters.html', 'w') as f:
    f.write(html)
