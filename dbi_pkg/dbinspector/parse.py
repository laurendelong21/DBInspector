import os
import os.path as osp
from dbinspector.startup import DATA, REFSEQ, REFSEQ_FASTA, UNIPROT
from dbinspector.utils import read_fasta, get_ncbi, get_uniprot
from time import time
from tqdm import tqdm
import logging
from typing import Dict

import gzip
import json
from collections import defaultdict
from lxml import etree
import pandas as pd


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def parse_all() -> None:
    """
    Wrapper function for parsing downloaded RefSeq and UniProt data that calls downstream functions, logs, tracks time.
    Results stored as json files in respective cache folders.
    """
    # ensure downloads are available
    download_data()
    # UniProt:
    parse_uniprot()
    # RefSeq in 3 steps:
    t0 = time()
    # (1) map RefSeq ID -> UniProt ID
    refseq_to_uniprot = map_refseq_to_uniprot()
    # (2) map RefSeq ID -> gene symbol
    refseq_to_symbol = map_refseq_to_symbol()
    # (3) assemble with sequences
    parse_refseq(refseq_to_uniprot, refseq_to_symbol)
    totaltime = (time() - t0)
    logger.info(f'...Finished parsing the RefSeq DB for human proteins in {totaltime:.2f} seconds.')
    logger.info("Parsing complete.")


def download_data():
    """
    Downloads database data if not already there.
    """
    # download necessary data via FTP
    if not (osp.exists(osp.join(DATA, "LRG_RefSeqGene"))
            and osp.exists(osp.join(DATA, "gene_refseq_uniprotkb_collab.gz"))):
        message = "Downloading the RefSeq data... this may take a few minutes."
        logger.info(message)
        print(message)
        get_ncbi('https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_refseq_uniprotkb_collab.gz', DATA)
        get_ncbi('https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene', DATA)
        for i in range(1, 9):
            get_ncbi(f'https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.{i}.protein.faa.gz', REFSEQ_FASTA)
    if not osp.exists(osp.join(DATA, "uniprot_sprot_human.xml.gz")):
        message = "Downloading the UniProt data... this may take a few minutes."
        logger.info(message)
        print(message)
        filename = 'uniprot_sprot_human.xml.gz'
        get_uniprot('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/'
                    + filename, DATA)


def map_refseq_to_uniprot() -> Dict[str, str]:
    """
    Reads the RefSeq_UniProt_collab file and maps RefSeq IDs to UniProt IDs.
    :return: dictionary of refseq: uniprot
    """
    ref_up_mapping = {}
    with gzip.open(osp.join(DATA, 'gene_refseq_uniprotkb_collab.gz'), 'rt') as filehandle:
        next(filehandle)  # skip header
        for line in tqdm(filehandle, desc='map RefSeq -> UniProt', leave=False):
            ref, up = line.split()
            # read in only relevant identifiers (approved proteins)
            if "NP_" in ref:
                ref_up_mapping[ref] = up
    return ref_up_mapping


def map_refseq_to_symbol() -> Dict[str, str]:
    """
    Gets the gene information for the RefSeq proteins based on the information available.
    :return: a dictionary of protein accessions.version : gene symbol
    """
    gentab = pd.read_csv(os.path.join(DATA, 'LRG_RefSeqGene'), sep='\t', usecols=['Protein', 'Symbol'])
    prot2symbol = {row['Protein']: [row['Symbol']] for i, row in gentab.iterrows()}
    return prot2symbol


def parse_uniprot() -> dict:
    """
    Process UniProt download into one dictionary as a return and to be saved as json.
    :return: dictionary of uniprot_id: {entry info}
    """
    t0, count = time(), 0
    namespace = '{http://uniprot.org/uniprot}'
    data = defaultdict(lambda: {'symbol': [],
                                'RefSeq ID': [],
                                'sequence': None})
    logger.info("Begin processing the human UniProt data ...")

    with gzip.open(osp.join(DATA, 'uniprot_sprot_human.xml.gz'), 'rb') as f:
        for event, elem in tqdm(etree.iterparse(f, events=("start", "end")), desc='parsing UniProt', leave=False):
            if elem.tag == namespace + 'accession':
                acc = elem.text
            if elem.tag == namespace + 'gene':
                for record in elem.findall(namespace + 'name'):
                    if record.get('type') == 'primary' or record.get('type') == 'synonym':
                        if record.text not in data[acc]['symbol']:
                            data[acc]['symbol'].append(record.text)
            if elem.tag == namespace + 'dbReference':
                if elem.attrib['type'] == 'RefSeq' and elem.attrib['id'] not in data[acc]['RefSeq ID']:
                    data[acc]['RefSeq ID'].append(elem.attrib['id'])
            if elem.tag == namespace + 'sequence':
                data[acc]['sequence'] = elem.text
            # delete parts of the tree to save memory
            while elem.getprevious() is not None:
                del elem.getparent()[0]  # clean up preceding siblings

    with open(osp.join(UNIPROT, f'uniprot.json'), 'w') as filehandle:
        json.dump(data, filehandle)
    totaltime = (time() - t0)
    logger.info(f'...Finished parsing the UniProt DB for human proteins in {totaltime:.2f} seconds.')
    return data


def parse_refseq(refseq_to_uniprot: Dict[str, str], refseq_to_symbol: Dict[str, str]) -> dict:
    """
    Final step in processing  RefSeq download into one dictionary as a return and to be saved as json
    :param dict refseq_to_uniprot: dictionary of refseq: uniprot
    :param dict refseq_to_symbol: a dictionary of protein accessions.version : gene symbol
    :return: dictionary of uniprot_id: {entry info}
    """
    logger.info("Attempting to begin processing the human RefSeq data...")
    data = defaultdict(lambda: {'symbol': [],
                                'UniProt ID': None,
                                'sequence': None})
    # process RefSeq sequences fasta file by fasta file
    for count, filename in enumerate(tqdm(os.listdir(REFSEQ_FASTA), desc='parsing RefSeq', leave=False)):
        if '.gz' not in filename:
            continue
        refseq_to_seq = read_fasta(osp.join(REFSEQ_FASTA, filename))  # (3) map RefSeq ID -> sequence

        # process the entries one by one
        for rsid, seq in refseq_to_seq.items():
            # RefSeq ID -> corresp. UniProt ID (1), symbol (2), sequence (3)
            data[rsid]['sequence'] = seq
            if rsid in refseq_to_symbol:
                data[rsid]['symbol'] = refseq_to_symbol[rsid]
            if rsid.split('.')[0] in refseq_to_uniprot:
                data[rsid]['UniProt ID'] = refseq_to_uniprot[rsid.split('.')[0]]

    with open(os.path.join(REFSEQ, f'refseq.json'), 'w') as filehandle:
        json.dump(data, filehandle)
    return data


if __name__ == '__main__':
    parse_all()
