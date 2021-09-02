import pytest
from dbinspector.utils import determine_identifier_type, format_list_entry, check_data_age, clear_dir
from dbinspector.utils import get_ncbi, get_uniprot
from dbinspector.exceptions import FileMissingError
import os
import os.path as osp
from pathlib import Path


HOME = str(Path.home())
TEST = osp.join(HOME, '.testfiles')
TEST_SUB = osp.join(TEST, 'testsubdir')
os.makedirs(TEST_SUB, exist_ok=True)


class TestUtils:
    """Test class for the frontend utils module"""
    def test_determine_identifer(self):
        """Test different identifier caching"""
        uniprot_id = determine_identifier_type('Q9NY95')
        refseq_id = determine_identifier_type('NP_149988.1')
        gene_symbol = determine_identifier_type('MICAR')
        empty = determine_identifier_type('')
        assert uniprot_id['uniprot_id']
        assert refseq_id['refseq_id']
        assert gene_symbol['symbol']
        assert not empty

    def test_format_list_entry(self):
        """Test if list is correctly formated"""
        str_list = ['a', 'b', 'c']
        assert format_list_entry(str_list) == 'a, b, c'
        assert format_list_entry([]) == ''

    def test_check_data_age(self):
        """Test if file age works properly"""
        with open(osp.join(TEST_SUB, 'test.txt'), 'w') as outfile:
            print("Test file", file=outfile)
        assert check_data_age(osp.join(TEST_SUB, 'test.txt'), "test file") == 0
        assert isinstance(check_data_age(osp.join(TEST_SUB, 'test.txt'), "test file"), int)
        with pytest.raises(FileMissingError):
            check_data_age(osp.join(TEST_SUB, 'test2.txt'), "2nd test file")

    def test_get_ncbi(self):
        """Test if file from NCBI FTP is successfully downloaded"""
        get_ncbi('https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/README', TEST)
        assert 'README' in set(os.listdir(TEST))
        # this test could be compromised if the README size is altered
        assert osp.getsize(osp.join(TEST, 'README')) > 6000

    def test_get_uniprot(self):
        """Test if file from UniProt FTP is successfully downloaded"""
        get_uniprot('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt', TEST)
        assert 'relnotes.txt' in set(os.listdir(TEST))
        # this test could be compromised if the relnotes.txt size is altered
        assert osp.getsize(osp.join(TEST, 'relnotes.txt')) > 900

    def test_clear_dir(self):
        """Test if the whole test directory is properly cleared without destroying the directory structure"""
        clear_dir(TEST)
        contents = os.listdir(TEST)
        dirs = [i for i in contents if osp.isdir(osp.join(TEST, i))]
        files = [i for i in contents if osp.isfile(osp.join(TEST, i))]
        assert len(dirs) == 1
        assert osp.join(TEST, dirs[0]) == TEST_SUB
        assert len(files) == 0
        assert len(os.listdir(TEST_SUB)) == 0
