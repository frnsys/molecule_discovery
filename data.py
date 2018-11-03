"""
Loaders for data
"""
import os
import re
import csv
import json
import sqlite3
import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
from lxml import etree
from scipy import sparse
from collections import defaultdict

here_ = os.path.join(os.path.dirname(__file__))
def here(p): return os.path.join(here_, p)

# Full schema:
# <ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_24_1_schema_documentation.txt>
CHEMBL_DB = here('data/chembl_24/chembl_24_sqlite/chembl_24.db')

def load_drugbank():
    tree = etree.parse(here('data/files/drugbank.xml'))
    ns = {'db': 'http://www.drugbank.ca'}

    drugs = []
    for drug in tqdm(tree.xpath('db:drug', namespaces=ns)):
        data = {}
        data['name'] = drug.find('db:name', ns).text

        data['ids'] = {}
        id_els = drug.findall('db:external-identifiers/db:external-identifier', ns)
        for el in id_els:
            source = el.find('db:resource', ns).text
            id = el.find('db:identifier', ns).text
            data['ids'][source] = id

        smiles = drug.find('db:calculated-properties/db:property[db:kind="SMILES"]/db:value', ns)
        if smiles is not None:
            smiles = smiles.text
        data['smiles'] = smiles

        codes = drug.findall('db:atc-codes/db:atc-code', ns)
        data['atcs'] = list(set([code.get('code') for code in codes]))

        data['targets'] = []
        target_els = drug.findall('db:targets/db:target', ns)
        for el in target_els:
            id = el.find('db:id', ns).text
            name = el.find('db:name', ns).text
            organism = el.find('db:organism', ns).text
            actions = [a.text for a in el.find('db:actions', ns)]
            uniprot_ids = [{
                'id': id.get('id'),
                'source': id.get('source')
            } for id in el.findall('db:polypeptide', ns)]

            data['targets'].append({
                'id': id,
                'uniprot_ids': uniprot_ids,
                'name': name,
                'organism': organism,
                'actions': actions,
            })
        drugs.append(data)

    # For debugging, seeing all target organisms
    # organisms = set()
    # for d in drugs:
    #     for t in d['targets']:
    #         organisms.add(t['organism'])
    # for o in organisms:
    #     print(o)

    # Filter to drugs with
    # human targets and CIDs and SMILES
    all_drugs = drugs
    drugs = []
    for d in all_drugs:
        cid = d['ids'].get('PubChem Compound')
        chembl = d['ids'].get('ChEMBL')
        if cid is None:
            continue

        # We filter to drugs that target humans;
        # note that this _does not_ include drugs that target
        # organisms/viruses that affect humans (e.g. HIV);
        # at some point we may want to include these as well
        targets = [t for t in d['targets']
                if t['organism'] is not None
                and t['organism'].lower() in ['human', 'homo sapiens']]
        if not targets:
            continue

        if d['smiles'] is None:
            continue

        drugs.append({
            'cid': cid,
            'chembl_id': chembl,
            'smiles': d['smiles'],
            'atcs': d['atcs'],
            'targets': targets
        })
    return drugs


def load_chembl():
    conn = sqlite3.connect(CHEMBL_DB)
    df = pd.read_sql_query('''
                        select
                            DRUG_MECHANISM.MECHANISM_OF_ACTION,
                            lower(DRUG_MECHANISM.ACTION_TYPE) as action,
                            COMPOUND_STRUCTURES.CANONICAL_SMILES as smiles,
                            MOLECULE_DICTIONARY.CHEMBL_ID as chembl_id,
                            COMPONENT_SEQUENCES.ACCESSION as uniprot_id,
                            COMPONENT_SEQUENCES.DB_SOURCE as source
                        from DRUG_MECHANISM
                        inner join COMPOUND_STRUCTURES
                            on DRUG_MECHANISM.MOLREGNO = COMPOUND_STRUCTURES.MOLREGNO
                        inner join MOLECULE_DICTIONARY
                            on DRUG_MECHANISM.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO
                        inner join TARGET_COMPONENTS
                            on DRUG_MECHANISM.TID = TARGET_COMPONENTS.TID
                        inner join COMPONENT_SEQUENCES
                            on TARGET_COMPONENTS.COMPONENT_ID = COMPONENT_SEQUENCES.COMPONENT_ID;
                        ''', conn)
    conn.close()
    return df.to_dict(orient='records')


def load_bindingdb():
    organism_col = 'Target Source Organism According to Curator or DataSource'
    df = pd.read_csv(here('data/bindingdb/BindingDB_All.tsv'), delimiter='\t', quoting=csv.QUOTE_NONE, error_bad_lines=False, dtype={'PubChem CID': 'str'})
    df = df[df[organism_col].str.contains('human|homo sapiens', case=False, na=False)]

    # Remove spaces from column names
    df.columns = df.columns.str.replace(r'\s+', '_')
    df.columns = df.columns.str.replace('.', '_')
    df.columns = df.columns.str.replace('(', '')
    df.columns = df.columns.str.replace(')', '')

    # Get UniProt target columns
    uniprot_col = 'UniProt_{}_Primary_ID_of_Target_Chain'
    swissprot_cols = 'SwissProt', [c for c in df.columns if c.startswith(uniprot_col.format('SwissProt'))]
    trembl_cols = 'TrEMBL', [c for c in df.columns if c.startswith(uniprot_col.format('TrEMBL'))]

    data = []
    for row in tqdm(df.itertuples()):
        uniprot_ids = []
        for src, cols in [swissprot_cols, trembl_cols]:
            for c in cols:
                id = getattr(row, c)
                if isinstance(id, str):
                    uniprot_ids.append({
                        'id': id,
                        'source': src,
                    })
        if not uniprot_ids:
            continue
        data.append({
            'cid': row.PubChem_CID,
            'ki': row.Ki_nM,
            'ec50': row.EC50_nM,
            'ic50': row.IC50_nM,
            'chembl_id': row.ChEMBL_ID_of_Ligand if isinstance(row.ChEMBL_ID_of_Ligand, str) else None,
            'targets': uniprot_ids,
            'smiles': row.Ligand_SMILES
        })
    return data


def load_compound_names():
    conn = sqlite3.connect(CHEMBL_DB)
    df = pd.read_sql_query('''
                        select
                            MOLECULE_DICTIONARY.CHEMBL_ID as chembl_id,
                            MOLECULE_DICTIONARY.PREF_NAME as name
                        from MOLECULE_DICTIONARY
                        where
                            MOLECULE_DICTIONARY.PREF_NAME is not null;
                        ''', conn)
    conn.close()
    chembl2cid = load_chembl2cid()
    return {chembl2cid.get(r['chembl_id'], r['chembl_id']): r['name'] for r in  df.to_dict(orient='records')}


def load_atc_codes():
    conn = sqlite3.connect(CHEMBL_DB)
    df = pd.read_sql_query('''
                        select
                            MOLECULE_DICTIONARY.CHEMBL_ID as chembl_id,
                            MOLECULE_ATC_CLASSIFICATION.LEVEL5 as atc_code
                        from MOLECULE_DICTIONARY
                        inner join MOLECULE_ATC_CLASSIFICATION
                            on MOLECULE_ATC_CLASSIFICATION.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO;
                        ''', conn)
    conn.close()
    chembl2cid = load_chembl2cid()
    atcs = defaultdict(set)
    for r in df.to_dict(orient='records'):
        id = chembl2cid.get(r['chembl_id'], r['chembl_id'])
        atcs[id].add(r['atc_code'])
    pubchem = json.load(open(here('data/files/pubchem_atc.json'), 'r'))
    for id, codes in pubchem.items():
        for code in codes:
            atcs[id].add(code)
    return atcs


def load_protein_names():
    df = pd.read_csv(here('data/files/uniprot_human.tsv'), delimiter='\t')
    df = df[['Entry', 'Protein names']]
    return {r['Entry']: r['Protein names'] for r in  df.to_dict(orient='records')}


def load_chembl2cid():
    chembl_re = re.compile('CHEMBL-?[0-9]+')
    lookup = {}
    with open(here('data/files/CID-CHEMBL'), 'r') as f:
        for line in tqdm(f):
            cid, part = line.strip().split(None, 1)
            chembl_id = chembl_re.search(part).group(0).replace('-', '')
            lookup[chembl_id] = cid
    return lookup


class Similarity:
    """Class to simplify similarity matrix lookups"""
    def __init__(self, mat, idx):
        self.mat = mat
        self.idx = idx

    def __getitem__(self, ids):
        a, b = ids
        try:
            i, j = sorted([self.idx[a], self.idx[b]])
            return self.mat[i,j]
        except KeyError:
            return None


def load_stitch():
    df = pd.read_csv(here('data/files/chemical_chemical.links.v5.0.tsv'), delimiter='\t')
    idx = {}

    # Extract CIDs
    for col in ['chemical1', 'chemical2']:
        df[col] = df[col].str[4:]

    # Prepare indices and similarity matrix
    for i, cid in enumerate(df['chemical1'].unique()):
        # Remove leading zeros
        cid = str(int(cid))
        idx[cid] = i
    n_cids = len(idx)
    mat = sparse.lil_matrix((n_cids, n_cids), dtype=np.float16)

    # Build similarity matrix
    for r in tqdm(df.itertuples()):
        cid1 = str(int(r.chemical1))
        cid2 = str(int(r.chemical2))
        i, j = sorted([idx[cid1], idx[cid2]])
        mat[i,j] = r.textmining/1000

    return Similarity(mat, idx)


def load_smiles(cids):
    smiles = {}
    for fn in tqdm(glob(here('data/smiles/*.smi'))):
        with open(fn, 'r') as f:
            for line in f:
                cid, smi = line.strip().split()
                if cid not in cids:
                    continue
                smiles[cid] = smi

    conn = sqlite3.connect(CHEMBL_DB)
    df = pd.read_sql_query('''
                        select
                            MOLECULE_DICTIONARY.CHEMBL_ID as chembl_id,
                            COMPOUND_STRUCTURES.CANONICAL_SMILES as smiles
                        from MOLECULE_DICTIONARY
                        inner join COMPOUND_STRUCTURES
                            on MOLECULE_DICTIONARY.MOLREGNO = COMPOUND_STRUCTURES.MOLREGNO;
                        ''', conn)
    conn.close()
    chembl_smiles = {}
    for r in df.to_dict(orient='records'):
        chembl_smiles[r['chembl_id']] = r['smiles']
    smiles.update(chembl_smiles)
    return smiles
