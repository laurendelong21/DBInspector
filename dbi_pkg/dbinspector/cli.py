import click
from dbinspector.compare import compare_entries, summary_statistics
from dbinspector.parse import parse_all
from dbinspector.utils import clear_dir, determine_identifier_type, check_data_age
import logging
import os.path as osp
import pandas as pd
from dbinspector.startup import DATA, UNIPROT, REFSEQ

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


@click.group()
def cli():
    """Entry method for the CLI."""
    pass


@cli.command()
def parse():
    """Parse the downloaded database data."""
    parse_all()
    click.echo("UniProt and RefSeq downloaded data parsed. Consider clearing large downloads with clear-cache.")


@cli.command()
@click.option("-q", "--query", type=str,
              help="UniProt ID, RefSeq ID, or gene symbol to compare entries in RefSeq and Uniprot.")
@click.option("-o", "--outfile", type=str, default=None,
              help="The filepath to which the comparison table should be written as tsv file.")
def compare(query: str, outfile: str = None):
    """
    Searches databases for matches of given query, compares entries across databases and prints results.
    Query can be UniProt ID, RefSeq ID, or gene symbol. Optionally saves results to a file
    :param str query: Database accession identifier or symbol to be compared.
    :param str outfile: Filepath for saving the results as a tsv.
    """
    if not query:
        params = {"uniprot_id": None, "refseq_id": None, "symbol": None}
    else:
        params = determine_identifier_type(query)
    comparison_df = compare_entries(refseq_id=params["refseq_id"],
                                    uniprot_id=params["uniprot_id"],
                                    symbol=params["symbol"])
    comparison_df.loc["sequence"] = comparison_df.loc["sequence"].str.wrap(30)
    pd.set_option('display.max_colwidth', 60, "display.max_columns", None, "expand_frame_repr", False)
    print(comparison_df)
    if outfile:
        comparison_df.to_csv(outfile, sep='\t')
        logger.info(f"comparison for {query} saved at {outfile}")

@cli.command()
@click.option("-o", "--outfile", type=str, default=None,
              help="The filepath to which the database summary should be written as tsv file, if desired.")
def database_summary(outfile: str = None):
    """
    Calculates summary statistics comparing the entries in the UniProt and RefSeq databases.
    Optionally saves results to a file.
    """
    stats_tab = summary_statistics()
    pd.set_option('display.max_colwidth', 60, "display.max_columns", None, "expand_frame_repr", False)
    print(stats_tab)
    if outfile:
        stats_tab.to_csv(outfile, sep='\t')
        logger.info(f"Database summary statistics saved at {outfile}")


@cli.command()
@click.option('-p', '--parsed', default=False, is_flag=True,
              help="A flag indicating that user wants to see age of parsed data.")
@click.option('-r', '--raw', default=False, is_flag=True,
              help="A flag indicating that user wants to see age of raw downloaded data.")
def check_age(parsed: bool = False, raw: bool = False):
    """
    Use this to check the age of raw and/or parsed files- if no flag specified, will show age of raw files.
    """
    if raw or (not raw and not parsed):
        check_data_age(osp.join(DATA, 'LRG_RefSeqGene'), 'raw refseq')
        check_data_age(osp.join(DATA, 'uniprot_sprot_human.xml.gz'), 'raw uniprot')
    if parsed:
        check_data_age(osp.join(REFSEQ, 'refseq.json'), 'parsed refseq')
        check_data_age(osp.join(UNIPROT, 'uniprot.json'), 'parsed uniprot')


@cli.command()
@click.option('-a', '--all_data', default=False, is_flag=True,
              help="Clear entire cache: Downloaded data as well as parsed files. Else only downloaded data.")
def clear_cache(all_data=False):
    """
    Use this to clear up the space used up by this program.
    Downloaded database data will be deleted, and if optional flag -a is used, so will the parsed files.
    """
    if all_data:
        if click.confirm('WARNING: Are you sure you want to delete all processed files?'
                         + ' Doing so will require parsing again to regain all functionality of this package.',
                         abort=True):
            click.echo('Clearing parsed data files from cache...')
            clear_dir(REFSEQ)
            clear_dir(UNIPROT)
    else:
        click.confirm('WARNING: Are you sure you want to delete the downloaded data files? '
                      + 'Project functionality remains, running parse after this will download/parse new, updated data',
                      abort=True)
    click.echo('Clearing downloaded data files from cache...')
    clear_dir(DATA)


if __name__ == '__main__':
    cli()
