import os.path as osp
from dbinspector.startup import REFSEQ, UNIPROT
import logging
from typing import Optional, Dict, List

import json

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# =======================================
#   for task 3.1: comprehensive mapping
# =======================================


def find_entries(refseq_id: str = None, uniprot_id: str = None,
                 symbol: str = None) -> Dict[str, List[Optional[dict]]]:
    """
    Wrapper function for mapping: Retrieves all available data for
     a given RefSeq ID, UniProt ID, or symbol.
    :return: all available data from each database
             {'RefSeq': dict, 'UniProt': dict}
             where each dict contains the corresponding ID of the
             other database, symbols, and amino acid sequences
    """
    if refseq_id:
        return retrieve_by_refseq_id(refseq_id)
    elif uniprot_id:
        return retrieve_by_uniprot_id(uniprot_id)
    elif symbol:
        return retrieve_by_symbol(symbol)


def retrieve_by_refseq_id(query: str) -> Dict[str, List[Optional[dict]]]:
    """
    Fetches all available information from cache on the given RefSeq
     identifier, including its corresponding UniProt entry.
    :return type: dict -- {'RefSeq': dict, 'UniProt': dict}
    """
    refseq_entry = get_refseq_entry(query)
    uniprot_entry = None
    if refseq_entry and refseq_entry['UniProt ID']:
        uniprot_entry = get_uniprot_entry(refseq_entry['UniProt ID'])
    return {'RefSeq': [refseq_entry] if refseq_entry else [],
            'UniProt': [uniprot_entry] if uniprot_entry else []}


def get_refseq_entry(query: str) -> Optional[dict]:
    """Retrieves a single RefSeq entry by accession ID."""
    data = read_refseq_data()
    try:
        refseq_entry = data[query]
        # give back the key as well
        refseq_entry['RefSeq ID'] = [query]
    except KeyError:
        refseq_entry = None
    return refseq_entry


def read_refseq_data() -> Dict[str, dict]:
    """
    Reads parsed RefSeq data from cache. Data is stored as a json file
     and read/returned as a dictionary.
    """
    with open(osp.join(REFSEQ, 'refseq.json')) as filehandle:
        data = json.load(filehandle)
    return data


def retrieve_by_uniprot_id(query: str) -> Dict[str, List[Optional[dict]]]:
    """
    Fetches all available information from cache on the given UniProt
     identifier, including corresponding RefSeq entries.
    :return type: dict -- {'UniProt': dict, 'RefSeq': list(dict)}
    """
    uniprot_entry = get_uniprot_entry(query)
    # for one UniProt entry there can be several corresponding RefSeq entries.
    refseq_entries = []
    if uniprot_entry and uniprot_entry['RefSeq ID']:
        for refseq_id in uniprot_entry['RefSeq ID']:
            refseq_entry = get_refseq_entry(refseq_id)
            if refseq_entry:
                refseq_entries.append(refseq_entry)
    return {'UniProt': [uniprot_entry] if uniprot_entry else [],
            'RefSeq': refseq_entries}


def get_uniprot_entry(query: str) -> Optional[dict]:
    """Retrieves a single UniProt entry by UniProt accession ID."""
    data = read_uniprot_data()
    try:
        uniprot_entry = data[query]
        # give back the key as well
        uniprot_entry['UniProt ID'] = query
    except KeyError:
        uniprot_entry = None
    return uniprot_entry


def read_uniprot_data() -> Dict[str, dict]:
    """Reads parsed UniProt data from cached json file."""
    with open(osp.join(UNIPROT, 'uniprot.json')) as filehandle:
        data = json.load(filehandle)
    return data


def retrieve_by_symbol(query: str) -> Dict[str, List[Dict[str, dict]]]:
    """
    Fetches all available entries from cached UniProt and RefSeq data
     on the given gene symbol.
    :return: '{UniProt': [{accession ID (str): {'symbol': (list(str)),
                                                'RefSeq ID': (list(str)),
                                                'sequence': (str)}}],
              'RefSeq':  [{accession ID (str):{'symbol': (list(str)),
                                               'UniProt ID': (str),
                                               'sequence': (str)}}] }
    """
    # UniProt
    uniprot_data, uniprot_matches = read_uniprot_data(), []
    for acc_id, entry in uniprot_data.items():
        if (query in entry['symbol']) or (query.upper() in entry['symbol']):
            match = entry.copy()
            # give back the key as well
            match['UniProt ID'] = acc_id
            uniprot_matches.append(match)
    # RefSeq
    refseq_data, refseq_matches = read_refseq_data(), []
    for acc_id, entry in refseq_data.items():
        if (query in entry['symbol']) or (query.upper() in entry['symbol']):
            match = entry.copy()
            # give back the key as well
            match['RefSeq ID'] = acc_id
            refseq_matches.append(match)
    return {'UniProt': uniprot_matches, 'RefSeq': refseq_matches}


if __name__ == '__main__':

    print("search by uniprot id:")
    print(find_entries(uniprot_id='Q9NSV2'))

    print("\nsearch by refseq id:")
    print(find_entries(refseq_id='NP_000585.2'))

    print("\nsearch by symbol:")
    temp = find_entries(symbol='ACTB')
    print('RefSeq', temp['RefSeq'], 'UniProt', temp['UniProt'], sep='\n')
