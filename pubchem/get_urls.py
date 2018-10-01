"""
Extract SDF download urls
from pubchem FTP directory
"""

import requests
import lxml.html
resp = requests.get('http://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/')
html = lxml.html.fromstring(resp.content)
urls = [a.attrib['href'] for a in html.cssselect('a')]
print(len(urls), 'urls')
with open('pubchem.txt', 'w') as f:
    f.write('\n'.join(urls))