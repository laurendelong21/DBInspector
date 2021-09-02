import pytest
import os.path as osp
from time import time
import json
from dbinspector.startup import DATA, REFSEQ, REFSEQ_FASTA, UNIPROT
from dbinspector.utils import clear_dir
from dbinspector.parse import download_data, map_refseq_to_uniprot, map_refseq_to_symbol
from dbinspector.parse import parse_uniprot, parse_refseq, parse_all


class TestParse:
    """Class for testing the parse commands."""
    def test_download_data(self):
        """Tests to see if the data downloads properly"""
        download_data()
        assert osp.exists(osp.join(DATA, "gene_refseq_uniprotkb_collab.gz"))
        assert osp.exists(osp.join(DATA, "LRG_RefSeqGene"))
        assert osp.exists(osp.join(REFSEQ_FASTA, "human.4.protein.faa.gz"))
        assert osp.exists(osp.join(DATA, "uniprot_sprot_human.xml.gz"))

    def test_map_refseq_to_uniprot(self):
        """Tests that mapping refseq IDs to uniprot IDs works"""
        refseq_to_uniprot = map_refseq_to_uniprot()
        assert refseq_to_uniprot['NP_001303897'] == 'Q96GV9'

    def test_map_refseq_to_symbol(self):
        """Tests that mapping refseq IDs to gene symbols works"""
        refseq_to_symbol = map_refseq_to_symbol()
        assert refseq_to_symbol['NP_000585.2'] == ['TNF']

    def test_parse_uniprot(self):
        """Tests that parsing the uniprot download into json metadata file works"""
        parse_uniprot()
        assert osp.exists(osp.join(UNIPROT, "uniprot.json"))
        assert osp.getsize(osp.join(UNIPROT, "uniprot.json")) > 13000000
        with open(osp.join(UNIPROT, "uniprot.json"), 'r') as uniprot_json:
            uniprot_dict = json.load(uniprot_json)
        assert uniprot_dict['Q9BWD0'] == {"symbol": ["ZWINT"],
                                          "RefSeq ID": ["NP_001005413.1", "NP_008988.2", "NP_127490.1"],
                                          "sequence": "MEAAETEAEAAALEVLAEVAGILEPVGLQEEAELPAKILVEFVVDSQKKDKLLCSQLQVADFLQ"
                                                      "NILAQEDTAKGLDPLASEDTSRQKAIAAKEQWKELKATYREHVEAIKIGLTKALTQMEEAQRKR"
                                                      "TQLREAFEQLQAKKQMAMEKRRAVQNQWQLQQEKHLQHLAEVSAEVRERKTGTQQELDRVFQKL"
                                                      "GNLKQQAEQERDKLQRYQTFLQLLYTLQGKLLFPEAEAEAENLPDDKPQQPTRPQEQSTGDTMG"
                                                      "RDPGVSFKAVGLQPAGDVNLP"}
        clear_dir(UNIPROT)

    def test_parse_refseq(self):
        """Tests that parsing the refseq download into json metadata file works"""
        t0 = time()
        refseq_to_uniprot = map_refseq_to_uniprot()  # (1) map RefSeq ID -> UniProt ID
        refseq_to_symbol = map_refseq_to_symbol()  # (2) map RefSeq ID -> gene symbol
        parse_refseq(refseq_to_uniprot, refseq_to_symbol)
        assert osp.exists(osp.join(REFSEQ, "refseq.json"))
        assert osp.getsize(osp.join(REFSEQ, "refseq.json")) > 39000000
        with open(osp.join(REFSEQ, "refseq.json"), 'r') as refseq_json:
            refseq_dict = json.load(refseq_json)
        assert refseq_dict['NP_001009958.1'] == {"symbol": ["ZNF655"],
                                                 "UniProt ID": "Q8N720",
                                                 "sequence": "MEEIPAQEAAGSPRVQFQSLETQSECLSPEPQFVQDTDMEQGLTGGILLRLPTTRI"
                                                             "HSVNSCPALSHTQASAFSGETLAVLTAGISKRWPKYRLPIDIARPCSETPFPRL"}
        clear_dir(REFSEQ)

    def test_parse_all(self):
        """Tests if the whole parse_all() pipeline works"""
        parse_all()
        assert osp.exists(osp.join(UNIPROT, "uniprot.json"))
        assert osp.getsize(osp.join(UNIPROT, "uniprot.json")) > 13000000
        assert osp.exists(osp.join(REFSEQ, "refseq.json"))
        assert osp.getsize(osp.join(REFSEQ, "refseq.json")) > 39000000
        clear_dir(UNIPROT)
        clear_dir(REFSEQ)

