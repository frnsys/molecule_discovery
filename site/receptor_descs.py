import json
import requests
import lxml.html
from tqdm import tqdm

BASE_URL = 'https://pharos.nih.gov/idg/targets/'

seen = set()

results = {}
clusters = json.load(open('data/clusters.json'))
for clus in tqdm(clusters):
    label = clus['label']
    receptor_id, action = label.split('__')
    if receptor_id in seen:
        continue

    resp = requests.get(BASE_URL + receptor_id)
    if resp.status_code != 200:
        print('Error for:', receptor_id)
        continue
    html = lxml.html.fromstring(resp.content)
    links = html.cssselect('a')
    link = next(l for l in links if 'models.Disease' in l.attrib.get('href', ''))
    url = 'https://pharos.nih.gov' + link.attrib['href']
    resp = requests.get(url)
    data = resp.json()
    results[receptor_id] = data

    seen.add(receptor_id)

with open('data/receptor_descs.json', 'w') as f:
    json.dump(results, f)
