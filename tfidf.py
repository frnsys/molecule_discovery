import math
from tqdm import tqdm
from collections import defaultdict

def tf_idf(docs):
    idf = defaultdict(int)
    tf = defaultdict(lambda: defaultdict(int))

    # Compute document frequencies
    # and term frequencies
    print('Computing DFs and TFs...')
    tok_docs = []
    for doc_id, toks in tqdm(docs):
        # Term counts
        for tok in toks:
            tf[doc_id][tok] += 1

        # Document frequencies
        # and set term frequencies
        for tok in set(toks):
            idf[tok] += 1
            tf[doc_id][tok] /= len(toks)

    # Compute (smoothed) inverse document frequencies
    print('Computing IDFs...')
    N = len(tok_docs) + 1
    for tok, count in tqdm(idf.items()):
        idf[tok] = math.log(1+(N/count))

    return tf, idf
