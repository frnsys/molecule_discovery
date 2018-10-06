"""
Note: these aren't used by may be useful
"""
import requests

UA = 'Mozilla/5.0 (X11; Linux x86_64; rv:64.0) Gecko/20100101 Firefox/64.0'

def get_pubmed(pmid):
    """Get PubMed info for a given PMID"""
    resp = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi', params={
        'db': 'pubmed',
        'retmode': 'json',
        'id': pmid
    })
    meta = resp.json()
    resp = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi', params={
        'db': 'pubmed',
        'retmode': 'text',
        'rettype': 'abstract',
        'id': pmid
    })
    abst = resp.content
    meta = meta['result'][str(pmid)]
    title = meta['title']
    return title, meta, abst

def get_uses(mol_id):
    """
    Note on request limits:

    > - No more than 5 requests per second.
    > - No more than 400 requests per minute.
    > - No longer than 300 second running time per minute.

    <https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access>
    """

    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON/?'.format(mol_id)
    resp = requests.get(url, headers={'User-Agent': UA})
    data = resp.json()
    sections = data['Record']['Section']
    uses = None
    for s in sections:
        if s['TOCHeading'] == 'Drug and Medication Information':
            for s_ in s['Section']:
                if s_['TOCHeading'] == 'Therapeutic Uses':
                    uses = s_['Information']
                    break
    return uses
