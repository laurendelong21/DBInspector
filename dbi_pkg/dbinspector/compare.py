from dbinspector.map import find_entries, read_uniprot_data, read_refseq_data
import pandas as pd
from dbinspector.exceptions import QueryNotFoundError, InputError
import dbinspector.startup
import logging
from typing import Dict, Optional, List

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def extract_query(arguments: List[str]) -> Optional[str]:
    """
    Extracts the query from the list of input arguments and asserts: exactly one query of data type string.
    :param list arguments: list of arguments containing the possible query and None
    :return: query, if no input error was triggered, i.e. query is conform.
    :rtype: str
    :raises InputError: if no query / too many queries / query is not a string
    """
    query = [q for q in arguments if q]
    if not query:
        logger.error("No query entered")
        raise InputError("No query entered")
    if len(query) != 1:
        logger.error('Too many queries were given')
        raise InputError('Too many queries were given')
    [query] = query  # unpack
    if not isinstance(query, str):
        logger.error('query must be type str')
        raise InputError('query must be type str')
    return query


def compare_entries(refseq_id: str = None, uniprot_id: str = None, symbol: str = None) -> pd.DataFrame:
    """Finds equivalent entries in RefSeq and UniProt based on a single query as a RefSeqID, UniProtID, or gene symbol.
    Results structured in a dataframe.

    :param str refseq_id: RefSeq ID to search for in the databases
    :param str uniprot_id: UniProt ID to search for in the databases
    :param str symbol: Gene symbol to search for in the databases
    :return: table displaying side-by-side matching entries
    :rtype: pd.DataFrame
    :raises QueryNotFoundError: if query not found in either database
    :raises InputError: if more than one or no queries entered.
    """
    query = extract_query([refseq_id, uniprot_id, symbol])
    logger.info(f"Searching for entries matching query: {query}")
    entries = find_entries(refseq_id, uniprot_id, symbol)
    if not entries["UniProt"] and not entries["RefSeq"]:  # query in neither
        logger.error(f"Query {query} could not be found")
        raise QueryNotFoundError(f"Query {query} as cannot be found in cached database download")

    if len(entries["UniProt"]) > 0:
        col_uniprot = [f"UniProt entry {i+1}" for i in range(len(entries["UniProt"]))]
    else:
        col_uniprot = ["No UniProt entry"]
        entries["UniProt"] = [{k: [] for k in entries["RefSeq"][0].keys()}]
    if len(entries["RefSeq"]) > 0:
        col_refseq = [f"RefSeq entry {i+1}" for i in range(len(entries["RefSeq"]))]
    else:
        col_refseq = ["No RefSeq entry"]
        entries["RefSeq"] = [{k: [] for k in entries["UniProt"][0].keys()}]
    columns = col_uniprot+col_refseq
    flattened_entries = entries["UniProt"] + entries["RefSeq"]
    entry_d = {'symbol': [],
         'UniProt ID': [],
         'RefSeq ID': [],
         'sequence': [],
         'sequence matches': [None]*len(flattened_entries),
         'sequence length': []}
    for i in range(len(flattened_entries)):
        entry_d['symbol'].append(flattened_entries[i]['symbol'])
        entry_d['UniProt ID'].append(flattened_entries[i]['UniProt ID'])
        entry_d['RefSeq ID'].append(flattened_entries[i]['RefSeq ID'])
        sequence = flattened_entries[i]['sequence']
        # check for sequence matches
        if sequence in entry_d['sequence']:
            all_index_matches = [index for index, seq in enumerate(entry_d['sequence']) if seq == sequence]
            all_index_matches.append(i)
            for j in all_index_matches:
                entry_d['sequence matches'][j] = [columns[k] for k in all_index_matches if k != j]
        entry_d['sequence'].append(sequence)
        entry_d['sequence length'].append(len(sequence))
    df = pd.DataFrame.from_dict(entry_d, orient="index", columns=columns)
    # pd.set_option("display.max_rows", None, "display.max_columns", None)
    logger.info(f"Found entries for {query}")
    return df


def summary_statistics() -> pd.DataFrame:
    """
    Compares and summarizes matches between metadata of human protein entries across databases.

    Checks matching gene symbol, UniProt accession ID, RefSeq accession ID, amino acid sequence and sequence length
    and gives a percentage of matches for each. The results are returned as a pandas DatFrame type.

    :return: table summarizing matching entries between databases
    :rtype: pd.DataFrame
    """
    logger.info("Producing database summary statistics")
    # uniprot often has several refseq equivalents
    # -> need to keep track of only unique entry matches or the percentage will be >100%
    consensus = {key: {"matches": 0, "UniProt entry": set()}
                 for key in ["Symbol", "RefSeq ID", "UniProt ID", "Sequence", "Sequence length"]}
    uniprot_data = read_uniprot_data()
    refseq_data = read_refseq_data()

    # find UniProt IDs in RefSeq entries - that are also IDs of UniProt entries
    for refseq_id in refseq_data:
        uniprot_id = refseq_data[refseq_id]["UniProt ID"]
        if uniprot_id and uniprot_id in uniprot_data:
            consensus = update_stats(refseq_data[refseq_id], uniprot_data[uniprot_id],
                                     refseq_id, uniprot_id, "refseq", consensus)

    # check the other direction: RefSeq IDs in UniProt entries in case this db shows equivalents that RefSeq doesn't
    for uniprot_id in uniprot_data:
        counterpart_ids = uniprot_data[uniprot_id]["RefSeq ID"]
        for refseq_id in counterpart_ids:
            if refseq_id in refseq_data and uniprot_id != refseq_data[refseq_id]["UniProt ID"]:
                # this only happens in one case
                consensus = update_stats(refseq_data[refseq_id], uniprot_data[uniprot_id],
                                         refseq_id, uniprot_id, "uniprot", consensus)

    return finalize_stats(consensus, len(uniprot_data), len(refseq_data))


def update_stats(refseq_data: dict, uniprot_data: dict, rsid: str, upid: str, first_db_searched: str, consensus: dict
                 ) -> Dict[str, int]:
    """Helper function used by summary_statistics(), not to be called by user."""
    # symbol match- refseq only has one (or no) symbol listed, uniprot has several (includng synonyms)
    if refseq_data['symbol'] and refseq_data['symbol'][0] in uniprot_data['symbol']:
        consensus["Symbol"]["matches"] += 1
        consensus["Symbol"]["UniProt entry"].add(upid)

    if first_db_searched == "refseq":
        # only called when uniprot ids do match
        consensus["UniProt ID"]["matches"] += 1
        consensus["UniProt ID"]["UniProt entry"].add(upid)
        if rsid in uniprot_data["RefSeq ID"]:
            consensus["RefSeq ID"]["matches"] += 1
            consensus["RefSeq ID"]["UniProt entry"].add(upid)

    elif first_db_searched == "uniprot":
        # only called when refseq ids do match
        consensus["RefSeq ID"]["matches"] += 1
        consensus["RefSeq ID"]["UniProt entry"].add(upid)
        if refseq_data["UniProt ID"] and upid in refseq_data["UniProt ID"]:
            consensus["UniProt ID"]["matches"] += 1
            consensus["UniProt ID"]["UniProt entry"].add(upid)

    if refseq_data['sequence'] == uniprot_data['sequence']:
        consensus["Sequence"]["matches"] += 1
        consensus["Sequence"]["UniProt entry"].add(upid)
        consensus["Sequence length"]["matches"] += 1
        consensus["Sequence length"]["UniProt entry"].add(upid)
    elif refseq_data['sequence'] and uniprot_data['sequence'] \
            and (len(refseq_data['sequence']) == len(uniprot_data['sequence'])):
        consensus["Sequence length"]["matches"] += 1
        consensus["Sequence length"]["UniProt entry"].add(upid)

    return consensus


def finalize_stats(consensus: Dict[str, int], num_entries_in_uniprot: int, num_entries_in_refseq: int) -> pd.DataFrame:
    """Helper function used by summary_statistics(), not to be called by user."""
    stats = {category: [value["matches"],
                        f"{len(value['UniProt entry'])/num_entries_in_uniprot:.2%}",
                        f"{value['matches']/num_entries_in_refseq:.2%}"]
             for category, value in consensus.items()}
    stats_df = pd.DataFrame.from_dict(stats, orient="index", columns=["Number of matches",
                                                                      "Matching UniProt entries [%]",
                                                                      "Matching RefSeq entries [%]"])
    return stats_df


if __name__ == '__main__':
    print(f"Summary Stats\n{'='*80}")
    print(summary_statistics(), '\n\n\n')

    ref1 = 'NP_149988.1'
    print(f"QUERY: {ref1}\n{'='*80}")
    print(compare_entries(refseq_id=ref1), '\n\n\n')

    # uniprot id which contains several refseq equivalents
    up1 = 'Q9NY95'
    print(f"QUERY: {up1}\n{'='*80}")
    print(compare_entries(uniprot_id=up1), '\n\n\n')

    # uniprot id which contains no refseq equivalents
    up2 = 'Q9NSV2'
    print(f"QUERY: {up2}\n{'=' * 80}")
    print(compare_entries(uniprot_id=up2), '\n\n\n')

    nonsense = "NOTaUNIPROTid"
    print(f"QUERY: {nonsense}\n{'='*80}")
    print(compare_entries(uniprot_id=nonsense),'\n')


