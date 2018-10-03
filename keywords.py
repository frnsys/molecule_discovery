import json
from tqdm import tqdm
from textblob import TextBlob


def get_keywords(text):
    tb = TextBlob(text)
    keywords = [t for t, tag in tb.tags if tag in ['NN', 'NNS', 'NNP', 'NNPS']]
    keywords.extend(tb.noun_phrases)
    return keywords


with open('pubmed.dat', 'r') as f:
    for line in tqdm(f.readlines()):
        data = json.loads(line.strip())
        kws = get_keywords(data['title'])
        kws.extend(get_keywords(data['abstract']))
        # print(kws)