import json
from tqdm import tqdm
from glob import glob
import pubmed_parser as pp

pmids = set()
with open('CID-PMID', 'r') as f:
    for line in tqdm(f.readlines()):
        line = line.strip()
        try:
            cid, pmid, _ = line.split()
            pmids.add(pmid)
        except:
            print(line)

datas = []
for path in tqdm(glob('pubmed/**/*.nxml')):
    article = pp.parse_pubmed_xml(path)
    pmid = article['pmid']
    if pmid not in pmids:
        continue
    data = {
        'pmid': pmid,
        'title': article['full_title'],
        'abstract': article['abstract']
    }
    datas.append(json.dumps(data))

print(len(datas), 'found')
with open('pubmed.dat', 'w') as f:
    f.write('\n'.join(datas))