from urllib import request
import logging
import os.path as osp
import ftputil
import os
import re
from typing import Union, Dict
import time
from dbinspector.exceptions import FileMissingError
import gzip


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def get_uniprot(url: str, download_dir: str) -> None:
    """Uses HTML request to access UniProt FTP download page and save as file.
    :param str url: should be the URL where the FTP download is located
    :param str download_dir: the place where the file should be stored
    :return: None; alternatively downloads the data in specified location
    """
    filename = osp.basename(url)
    logger.info(f'Beginning download from {filename}...')
    request.urlretrieve(url, osp.join(download_dir, filename))
    logger.info(f'...Finished file download for {filename}.')


def get_ncbi(url: str, download_dir: str) -> None:
    """Uses HTML request to access NCBI FTP download page and save as file.
    :param str url: should be the URL where the FTP download is located, including the file name at the end
    :param str download_dir: the place where the file should be stored
    :return: None; alternatively downloads the data in specified location
    """
    url = url.replace('https://', '')
    hostdir = osp.dirname(url).replace('ftp.ncbi.nlm.nih.gov', '') + '/'
    filename = osp.basename(url)
    host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')
    host.chdir(hostdir)
    logger.info(f'Beginning download from {filename}...')
    host.download(filename, osp.join(download_dir, filename))
    logger.info(f'...Finished file download for {filename}.')


def read_fasta(filename: str) -> Dict[str, str]:
    """
    Reads a file containing multiple FASTA sequences and returns a dictionary of the header: sequence
    :param str filename: should be the name of the open file to be read
    :return: dict containing the header: the sequence
    """
    seq_dict = {}
    with gzip.open(filename, 'rt') as fc:
        all_lines = str(fc.read())
        seqs = all_lines.split('>')
        for seq in seqs[1:]:
            seq = seq.strip('"').split('\n')
            ref_id, prot_name = seq[0].replace(' [Homo sapiens]', '').split(' ', 1)
            if "NP_" in ref_id:
                seq_dict[ref_id] = ''.join(seq[1:])
        return seq_dict


def clear_dir(directory: str) -> None:
    """Recursively clears all files in given directory and its subdirectories
    :param str directory: the directory through which to recurse and delete all files
        (will maintain subdirectory structure but delete all files in subdirectories)
    """
    logger.info(f"Clearing all files in the {directory} directory...")
    for item in os.listdir(directory):
        path = osp.join(directory, item)
        if not osp.isdir(path):
            os.remove(path)
        else:
            clear_dir(path)


def determine_identifier_type(identifier: str) -> Union[dict, None]:
    """
    Can be used to determine the type of identifier.
    :param str identifier: String input of an refseq_id, uniprot_id or gene symbol
    :return: dictionary containing the most likely mapping between the identifier and the possible id type
    """
    uniprot_id_pattern = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    if not identifier:
        return None
    if "NP_" in identifier:
        return {"uniprot_id": None, "refseq_id": identifier, "symbol": None}
    elif re.match(uniprot_id_pattern, identifier):
        return {"uniprot_id": identifier, "refseq_id": None, "symbol": None}
    else:
        return {"uniprot_id": None, "refseq_id": None, "symbol": identifier}


def format_list_entry(str_list: list) -> str:
    """
    Converts list of strings into a single string with values separated by comma.
    :param list str_list: list of strings
    :return: formatted string
    """
    out_str = ''
    for elm in str_list:
        out_str += elm
        out_str += ', '
    return out_str[:-2]


def check_data_age(file_path: str, data_name: str) -> int:
    """Checks to see if the file is over 1 week old, and if so, asks the user if they want to refresh the data.
    :param str file_path: the filepath to the file to check
    :param str data_name: a name of the data
    :return: rounded age of file in days
    :raises FileMissingError: if file in question does not exist
    """
    if osp.exists(file_path):
        file_age_sec = time.time() - osp.getmtime(file_path)
        file_age_days = round(file_age_sec / 86400)
        logger.info(f"File age for {file_path} was determined to be {file_age_days} day(s) old.")
        print(f"Your {data_name} data is {file_age_days} day(s) old.")
        return file_age_days
    else:
        logger.warning(f"The age of {file_path} can not be determined because this file or directory does not exist.")
        raise FileMissingError
