"""
Extract SDF download urls
from pubchem FTP directory
"""
import requests
import lxml.html

def get_urls(base):
    resp = requests.get('http' + base)
    html = lxml.html.fromstring(resp.content)
    urls = [a.attrib['href'] for a in html.cssselect('a')]
    urls = ['ftp' + base + url for url in urls if url.endswith('.gz')]
    return urls

# Get PubChem SDF URLs
base = '://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/'
urls = get_urls(base)
print(len(urls), 'urls')
with open('pubchem.txt', 'w') as f:
    f.write('\n'.join(urls))

# Get PubMed Open Access URLs
base = '://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/'
urls = get_urls(base)
print(len(urls), 'urls')
with open('pubmed.txt', 'w') as f:
    f.write('\n'.join(urls))