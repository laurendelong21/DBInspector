import json
import pandas as pd
import os
import os.path as osp
from pathlib import Path

from dbinspector.map import read_refseq_data, read_uniprot_data, get_refseq_entry, get_uniprot_entry, \
    retrieve_by_refseq_id, retrieve_by_uniprot_id, retrieve_by_symbol, find_entries
from dbinspector.compare import summary_statistics

from dbinspector.startup import REFSEQ, UNIPROT, CACHE
REFSEQ_DATA = osp.join(REFSEQ, 'refseq.json')
REFSEQ_TEMP = osp.join(CACHE, 'refseq.json')
UNIPROT_DATA = osp.join(UNIPROT, 'uniprot.json')
UNIPROT_TEMP = osp.join(CACHE, 'uniprot.json')


def metadata_keys_complete(result_keys) -> bool:
    """
    Helper function for tests that assert, that all metdata keys are present in a dictionary that has been
    returned as the result of a query.

    :param result_keys: keys of a database-specific dictionary that has been returned by a query, corresponds to
                        metadata categories
    :type result_keys: dict_keys
    """
    return set(result_keys) == {'symbol', 'UniProt ID', 'sequence', 'RefSeq ID'}


class TestEnvironment:
    def test_environment(self):
        """Makes sure that the test data is available which is the underlying assumptions of functions in map."""
        # if available - move parsed data up
        if osp.exists(REFSEQ_DATA) and osp.exists(UNIPROT_DATA):
            Path(REFSEQ_DATA).rename(REFSEQ_TEMP)
            Path(UNIPROT_DATA).rename(UNIPROT_TEMP)
            assert osp.isfile(REFSEQ_TEMP)
            assert osp.getsize(REFSEQ_TEMP) > 39000000
            assert osp.isfile(UNIPROT_TEMP)
            assert osp.getsize(UNIPROT_TEMP) > 13000000
        # replace data by test data
        create_test_data()
        assert osp.isfile(REFSEQ_DATA)
        assert osp.getsize(REFSEQ_DATA) == 832
        assert osp.isfile(UNIPROT_DATA)
        assert osp.getsize(UNIPROT_DATA) == 702


class TestMap:
    """Tests the dbinspector modules map and compare on static test data."""

    def test_read_parsed_data(self):
        """Checks that parsed test data is loaded/read in correctly from the cached json files to a dictionary."""
        refseq = read_refseq_data()
        uniprot = read_uniprot_data()
        # data structure
        assert isinstance(refseq, dict)
        assert all([isinstance(key, str) for key in refseq.keys()])
        assert all([isinstance(val, dict) for val in refseq.values()])
        assert isinstance(uniprot, dict)
        assert all([isinstance(key, str) for key in uniprot.keys()])
        assert all([isinstance(val, dict) for val in uniprot.values()])
        # number of entries
        assert len(refseq) == 11
        assert len(uniprot) == 9

    def test_get_refseq_entry(self):
        """Tests the lookup of a RefSeq entry by accession ID by comparison to an example."""
        res: dict = get_refseq_entry('rsid1')
        assert metadata_keys_complete(res.keys())
        assert res['symbol'] == ['ONE']
        assert res['UniProt ID'] == 'upid1'
        assert res['sequence'] == 'ONEONEONE'
        assert res['RefSeq ID'] == ['rsid1']

    def test_get_uniprot_entry(self):
        """Tests the lookup of a UniProt entry by accession ID by comparison to an example."""
        res: dict = get_uniprot_entry('upid1')
        assert metadata_keys_complete(res.keys())
        assert res['symbol'] == ['ONE']
        assert res['UniProt ID'] == 'upid1'
        assert res['sequence'] == 'ONEONEONE'
        assert res['RefSeq ID'] == ['rsid1']

    def test_retrieve_by_refseq_id(self):
        """Tests retrieval of appropriate RefSeq and UniProt entries upon query by RefSeq accession ID."""
        # ==== successful ==== 1 RefSeq & 1 UniProt entry ================
        res: dict = retrieve_by_refseq_id('rsid2')
        assert res['RefSeq'] and res['UniProt']
        # RefSeq
        assert isinstance(res['RefSeq'], list)
        assert len(res['RefSeq']) == 1
        assert isinstance(res['RefSeq'][0], dict)
        assert metadata_keys_complete(res['RefSeq'][0].keys())
        assert res['RefSeq'][0]['symbol'] == ['TWO']
        assert res['RefSeq'][0]['UniProt ID'] == 'upid2'
        assert res['RefSeq'][0]['sequence'] == 'TWOTWOTWO'
        assert res['RefSeq'][0]['RefSeq ID'] == ['rsid2']
        # UniProt
        assert isinstance(res['UniProt'], list)
        assert len(res['UniProt']) == 1
        assert isinstance(res['UniProt'][0], dict)
        assert metadata_keys_complete(res['UniProt'][0].keys())
        assert res['UniProt'][0]['symbol'] == ['ZWEI']
        assert res['UniProt'][0]['UniProt ID'] == 'upid2'
        assert res['UniProt'][0]['sequence'] == 'ZWEIZWEIZWEI'
        assert res['UniProt'][0]['RefSeq ID'] == ['rsid2']
        # ==== semi-successful == 1 RefSeq & No UniProt entry ============
        res: dict = retrieve_by_refseq_id('rsid6')
        assert res['RefSeq'] and not res['UniProt']
        assert isinstance(res['RefSeq'], list)
        assert len(res['RefSeq']) == 1
        assert isinstance(res['RefSeq'][0], dict)
        assert metadata_keys_complete(res['RefSeq'][0].keys())
        assert res['RefSeq'][0]['symbol'] == ['SIX']
        assert not res['RefSeq'][0]['UniProt ID']
        assert res['RefSeq'][0]['sequence'] == 'SIXSIX'
        assert res['RefSeq'][0]['RefSeq ID'] == ['rsid6']

    def test_retrieve_by_uniprot_id(self):
        """Tests retrieval of appropriate RefSeq and UniProt entries upon query by RefSeq accession ID."""
        # ==== successful ==== 1 RefSeq & 1 UniProt entry ================
        res: dict = retrieve_by_uniprot_id('upid2')
        assert res['RefSeq'] and res['UniProt']
        # RefSeq
        assert isinstance(res['RefSeq'], list)
        assert len(res['RefSeq']) == 1
        assert isinstance(res['RefSeq'][0], dict)
        assert metadata_keys_complete(res['RefSeq'][0].keys())
        assert res['RefSeq'][0]['symbol'] == ['TWO']
        assert res['RefSeq'][0]['UniProt ID'] == 'upid2'
        assert res['RefSeq'][0]['sequence'] == 'TWOTWOTWO'
        assert res['RefSeq'][0]['RefSeq ID'] == ['rsid2']
        # UniProt
        assert isinstance(res['UniProt'], list)
        assert len(res['UniProt']) == 1
        assert isinstance(res['UniProt'][0], dict)
        assert metadata_keys_complete(res['UniProt'][0].keys())
        assert res['UniProt'][0]['symbol'] == ['ZWEI']
        assert res['UniProt'][0]['UniProt ID'] == 'upid2'
        assert res['UniProt'][0]['sequence'] == 'ZWEIZWEIZWEI'
        assert res['UniProt'][0]['RefSeq ID'] == ['rsid2']
        # ==== semi-successful == 1 UniProt & No RefSeq entry ============
        res: dict = retrieve_by_uniprot_id('upid6')
        assert res['UniProt'] and not res['RefSeq']
        assert isinstance(res['UniProt'], list)
        assert len(res['UniProt']) == 1
        assert isinstance(res['UniProt'][0], dict)
        assert metadata_keys_complete(res['UniProt'][0].keys())
        assert res['UniProt'][0]['symbol'] == ['SIX']
        assert res['UniProt'][0]['UniProt ID'] == 'upid6'
        assert res['UniProt'][0]['sequence'] == 'SIXSIX'
        assert not res['UniProt'][0]['RefSeq ID']
        # ==== successful ==== 1 UniProt & 2 RefSeq entries ==============
        res = retrieve_by_uniprot_id('upid4')
        assert res['RefSeq'] and res['UniProt']
        # UniProt
        assert isinstance(res['UniProt'], list)
        assert len(res['UniProt']) == 1
        assert isinstance(res['UniProt'][0], dict)
        assert metadata_keys_complete(res['UniProt'][0].keys())
        assert set(res['UniProt'][0]['symbol']) == {'FOUR', 'CUATRO'}
        assert res['UniProt'][0]['UniProt ID'] == 'upid4'
        assert res['UniProt'][0]['sequence'] == 'CUATRO44'
        assert set(res['UniProt'][0]['RefSeq ID']) == {'rsid4', 'rsid5'}
        # RefSeq
        assert isinstance(res['RefSeq'], list)
        assert len(res['RefSeq']) == 2
        #     entry rsid4
        assert isinstance(res['RefSeq'][0], dict)
        assert metadata_keys_complete(res['RefSeq'][0].keys())
        assert res['RefSeq'][0]['symbol'] == ['FOUR']
        assert res['RefSeq'][0]['UniProt ID'] == 'upid4'
        assert res['RefSeq'][0]['sequence'] == 'FOURFOUR'
        assert res['RefSeq'][0]['RefSeq ID'] == ['rsid4']
        #     entry rsid5
        assert isinstance(res['RefSeq'][1], dict)
        assert metadata_keys_complete(res['RefSeq'][1].keys())
        assert not res['RefSeq'][1]['symbol']
        assert not res['RefSeq'][1]['UniProt ID']
        assert res['RefSeq'][1]['sequence'] == 'FIVEFIVE'
        assert res['RefSeq'][1]['RefSeq ID'] == ['rsid5']

    def test_retrieve_by_symbol(self):
        """Tests that retrieval by symbol returns entries from UniProt and RefSeq in the appropriate data structure."""
        # ==== semi-successful == No UniProt & 1 RefSeq entry ============
        res: dict = retrieve_by_symbol('TWO')
        assert res['RefSeq'] and not res['UniProt']
        assert isinstance(res['RefSeq'], list)
        assert len(res['RefSeq']) == 1
        assert isinstance(res['RefSeq'][0], dict)
        assert metadata_keys_complete(res['RefSeq'][0].keys())
        assert res['RefSeq'][0]['symbol'] == ['TWO']
        assert res['RefSeq'][0]['UniProt ID'] == 'upid2'
        assert res['RefSeq'][0]['sequence'] == 'TWOTWOTWO'
        assert res['RefSeq'][0]['RefSeq ID'] == 'rsid2'
        # ==== successful ==== 2 UniProt & 1 RefSeq entries ==============
        res: dict = retrieve_by_symbol('EIGHT')
        assert res['RefSeq'] and res['UniProt']
        # UniProt
        assert isinstance(res['UniProt'], list)
        assert len(res['UniProt']) == 2
        #     entry upid7
        assert isinstance(res['UniProt'][0], dict)
        assert metadata_keys_complete(res['UniProt'][0].keys())
        assert res['UniProt'][0]['symbol'] == ['EIGHT']
        assert res['UniProt'][0]['UniProt ID'] == 'upid7'
        assert not res['UniProt'][0]['sequence']
        assert not res['UniProt'][0]['RefSeq ID']
        #     entry upid8
        assert isinstance(res['UniProt'][1], dict)
        assert metadata_keys_complete(res['UniProt'][1].keys())
        assert set(res['UniProt'][1]['symbol']) == {'EIGHT', 'OCHO'}
        assert res['UniProt'][1]['UniProt ID'] == 'upid8'
        assert not res['UniProt'][1]['sequence']
        assert not res['UniProt'][1]['RefSeq ID']
        # RefSeq
        assert isinstance(res['RefSeq'], list)
        assert len(res['RefSeq']) == 1
        assert isinstance(res['RefSeq'][0], dict)
        assert metadata_keys_complete(res['RefSeq'][0].keys())
        assert set(res['RefSeq'][0]['symbol']) == {'EIGHT', 'ACHT'}
        assert not res['RefSeq'][0]['UniProt ID']
        assert not res['RefSeq'][0]['sequence']
        assert res['RefSeq'][0]['RefSeq ID'] == 'rsid8'

    def test_find_entries(self):
        """Tests the appropriate functions' return data types."""
        assert isinstance(find_entries(refseq_id='rsid6'), dict)
        assert isinstance(find_entries(refseq_id='rsid6')['RefSeq'], list)
        assert isinstance(find_entries(uniprot_id='upid6'), dict)
        assert isinstance(find_entries(uniprot_id='upid6')['UniProt'], list)


class TestCompare:
    """Testing class to check proper functionality of methods in module dbinspector.compare on test data"""
    def test_summary_statistics(self) -> None:
        """Determines the correctness of calculated summary statistics on the generated test data."""
        res = summary_statistics()
        # data structure
        assert isinstance(res, pd.DataFrame)
        assert set(res.columns) == {'Number of matches', 'Matching UniProt entries [%]', 'Matching RefSeq entries [%]'}
        assert set(res.index) == {'Symbol', 'RefSeq ID', 'UniProt ID', 'Sequence', 'Sequence length'}
        # values
        matches = res['Number of matches']
        assert matches['Symbol'] == 4
        assert matches['RefSeq ID'] == 7
        assert matches['UniProt ID'] == 7
        assert matches['Sequence'] == 3
        assert matches['Sequence length'] == 5
        uniprot = res['Matching UniProt entries [%]']
        assert uniprot['Symbol'] == '44.44%'
        assert uniprot['RefSeq ID'] == '44.44%'
        assert uniprot['UniProt ID'] == '55.56%'
        assert uniprot['Sequence'] == '33.33%'
        assert uniprot['Sequence length'] == '44.44%'
        refseq = res['Matching RefSeq entries [%]']
        assert refseq['Symbol'] == '36.36%'
        assert refseq['RefSeq ID'] == '63.64%'
        assert refseq['UniProt ID'] == '63.64%'
        assert refseq['Sequence'] == '27.27%'
        assert refseq['Sequence length'] == '45.45%'


class TestEnvironmentRestore:
    def test_environment_exit(self):
        """Exits the test data mode by replacing the test data with parsed files again."""
        # remove test data
        os.remove(REFSEQ_DATA)
        os.remove(UNIPROT_DATA)
        # if available - move back parsed data
        if osp.exists(REFSEQ_TEMP) and osp.exists(UNIPROT_TEMP):
            Path(REFSEQ_TEMP).rename(REFSEQ_DATA)
            Path(UNIPROT_TEMP).rename(UNIPROT_DATA)
            assert osp.getsize(REFSEQ_DATA) > 39000000
            assert osp.getsize(UNIPROT_DATA) > 13000000
        assert not osp.isfile(REFSEQ_TEMP)
        assert not osp.isfile(UNIPROT_TEMP)


def create_test_data() -> None:
    refseq = {'rsid1': {'symbol': ['ONE'], 'UniProt ID': 'upid1', 'sequence': 'ONEONEONE'},
              'rsid2': {'symbol': ['TWO'], 'UniProt ID': 'upid2', 'sequence': 'TWOTWOTWO'},
              'rsid3': {'symbol': ['THREE'], 'UniProt ID': 'upid3', 'sequence': 'THREETHREE'},
              'rsid4': {'symbol': ['FOUR'], 'UniProt ID': 'upid4', 'sequence': 'FOURFOUR'},
              'rsid5': {'symbol': [], 'UniProt ID': None, 'sequence': 'FIVEFIVE'},
              'rsid6': {'symbol': ['SIX'], 'UniProt ID': None, 'sequence': 'SIXSIX'},
              'rsid7': {'symbol': [], 'UniProt ID': None, 'sequence': ''},
              'rsid8': {'symbol': ['EIGHT', 'ACHT'], 'UniProt ID': None, 'sequence': ''},
              'rsid9': {'symbol': ['NEUN'], 'UniProt ID': 'upid9', 'sequence': 'NEUNNEUN'},
              'rsid10': {'symbol': ['NEUF'], 'UniProt ID': 'upid9', 'sequence': 'NUEVENUEVE'},
              'rsid11': {'symbol': ['NUEVE'], 'UniProt ID': 'upid9', 'sequence': 'NINENINE'}}

    uniprot = {'upid1': {'symbol': ['ONE'], 'RefSeq ID': ['rsid1'], 'sequence': 'ONEONEONE'},
               'upid2': {'symbol': ['ZWEI'], 'RefSeq ID': ['rsid2'], 'sequence': 'ZWEIZWEIZWEI'},
               'upid3': {'symbol': ['THREE'], 'RefSeq ID': [], 'sequence': 'THREETHREE'},
               'upid4': {'symbol': ['FOUR', 'CUATRO'], 'RefSeq ID': ['rsid4', 'rsid5'], 'sequence': 'CUATRO44'},
               'upid5': {'symbol': [], 'RefSeq ID': [], 'sequence': ''},
               'upid6': {'symbol': ['SIX'], 'RefSeq ID': [], 'sequence': 'SIXSIX'},
               'upid7': {'symbol': ['EIGHT'], 'RefSeq ID': [], 'sequence': ''},
               'upid8': {'symbol': ['EIGHT', 'OCHO'], 'RefSeq ID': [], 'sequence': ''},
               'upid9': {'symbol': ['NUEVE'], 'RefSeq ID': ['rsid9', 'rsid10', 'rsid11'], 'sequence': 'NUEVENUEVE'}}
    # save
    with open(REFSEQ_DATA, 'w') as filehandle:
        json.dump(refseq, filehandle)
    with open(UNIPROT_DATA, 'w') as filehandle:
        json.dump(uniprot, filehandle)
