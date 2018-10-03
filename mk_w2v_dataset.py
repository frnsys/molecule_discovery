# python -m spacy download en_core_web_sm
import math
import spacy
from tqdm import tqdm
from glob import glob
import pubmed_parser as pp
from multiprocessing import Pool, cpu_count

n_cores = cpu_count()
chunk_size = 10
nlp = spacy.load('en_core_web_sm')
pm_files = glob('pubmed/**/*.nxml')

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def stream(files):
    for path in files:
        article = pp.parse_pubmed_xml(path)
        text = '. '.join([
            article['full_title'],
            article['abstract']
        ])
        yield text

def extract_tokens(files):
    tokens = []
    for text in stream(files):
        doc = nlp(text)
        tokens.append('<START>')
        tokens.extend(tok.text for tok in doc)
        tokens.append('<END>')
    return tokens

p = Pool(n_cores)
with open('output.txt', 'a') as f:
    total = math.ceil(len(pm_files)/chunk_size)
    for tokens in tqdm(p.imap(extract_tokens, chunks(pm_files, chunk_size)), total=total):
        f.write(' '.join(tokens))
        f.write(' ')
p.close()
