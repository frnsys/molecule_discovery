"""
Download ATC codes from PubChem
"""
import re
import json
import requests
from tqdm import tqdm
from collections import defaultdict

CODE_RE = re.compile('code=([A-Z0-9]+)')
URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/JSON/?source=WHO%20ATC&heading=ATC+Code&response_type=save&response_basename=PubChemAnnotations_source=WHO%20ATC&heading=ATC+Code&page={page}'

def extract_atc_code(d):
    try:
        cid = d['LinkToPubChemBy']['CID'][0]
    except KeyError:
        return None
    url = d['Data'][0]['URL']
    code = CODE_RE.search(url).group(1)
    return cid, code

raw_atcs = []
resp = requests.get(URL.format(page=1))
data = resp.json()['Annotations']
n_pages = data['TotalPages']
raw_atcs += data['Annotation']

for i in tqdm(range(2, n_pages+1)):
    resp = requests.get(URL.format(page=i))
    data = resp.json()['Annotations']
    raw_atcs += data['Annotation']

atcs = defaultdict(list)
for cid, code in filter(None, map(extract_atc_code, raw_atcs)):
    atcs[cid].append(code)

print('ATC annotations:', len(atcs))
with open('files/pubchem_atc.json', 'w') as f:
    json.dump(atcs, f)
