from dbinspector.parse import parse_all
from dbinspector.map import get_refseq_entry, get_uniprot_entry
from dbinspector.compare import extract_query, compare_entries, summary_statistics, update_stats, finalize_stats

from dbinspector.startup import REFSEQ, UNIPROT
import os.path as osp

import pandas as pd

NONSENSE = 'Blubb~Blubb~Blubb~O0o.'
MSG = "This test could be compromised because of updated metadata from the database."


def is_approx_percentage(result: str, comparison: float) -> bool():
    """Determines whether or not a value, given in % as a string, is less than 5% higher/lower than the comparison."""
    result = float(result[:-1])
    return 0.95 * comparison < result < 1.05 * comparison


def is_approx(result: str, comparison: int) -> bool():
    """Determines whether or not a value, given as string, is less than 5% higher/lower than the comparison."""
    return 0.95 * comparison < int(result) < 1.05 * comparison


class TestCompare:
    """Testing class to check proper functionality of methods in module dbinspector.compare"""

    def test_preparation(self):
        """Makes sure that the parsed data is available which is the underlying assumptions of functions in compare."""
        parse_all()
        refseq_path = osp.join(REFSEQ, 'refseq.json')
        uniprot_path = osp.join(UNIPROT, 'uniprot.json')
        assert osp.isfile(refseq_path)
        assert osp.isfile(uniprot_path)

    def test_extract_query(self) -> None:
        """Tests that for a correct query the assertion of the query type runs without assertion errors."""
        assert extract_query(['NP_149988.1', None, None]) == 'NP_149988.1'
        assert extract_query([None, 'Q9NSV2', None]) == 'Q9NSV2'
        assert extract_query([None, None, 'ACTB']) == 'ACTB'
        assert extract_query([NONSENSE])

    def test_compare_entries_refseqid(self) -> None:
        """Tests the comparison of database entries using a given example from a successful query with RefSeq ID."""
        res = compare_entries(refseq_id='NP_149988.1', uniprot_id=None, symbol=None)
        # data structure
        assert isinstance(res, pd.DataFrame)
        assert set(res.columns) == {'UniProt entry 1', 'RefSeq entry 1'}, MSG
        assert set(res.index) == {'symbol', 'UniProt ID', 'RefSeq ID', 'sequence', 'sequence matches',
                                  'sequence length'}
        # values UniProt
        seq = ('MEVDINGESRSTLTTLPFPGAEANSPGKAEAEKPRCSSTPCSPMRRTVSGYQILHMDSNYLVGFTTGEEL'
               + 'LKLAQKCTGGEESKAEAMPSLRSKQLDAGLARSSRLYKTRSRYYQPYEIPAVNGRRRRRMPSSGDKCT'
               + 'KSLPYEPYKALHGPLPLCLLKGKRAHSKSLDYLNLDKMIKEPADTEVLQYQLQHLTLRGDRVFARNNT')
        assert set(res['UniProt entry 1']['symbol']) == {'MACIR', 'C5orf30'}, MSG
        assert res['UniProt entry 1']['UniProt ID'] == 'Q96GV9', MSG
        assert set(res['UniProt entry 1']['RefSeq ID']) == {'NP_001303897.1', 'NP_001303898.1', 'NP_149988.1'}, MSG
        assert res['UniProt entry 1']['sequence'] == seq, MSG
        assert res['UniProt entry 1']['sequence matches'] == ['RefSeq entry 1'], MSG
        assert res['UniProt entry 1']['sequence length'] == 206, MSG
        # values RefSeq
        assert not res['RefSeq entry 1']['symbol'], MSG
        assert res['RefSeq entry 1']['UniProt ID'] == 'Q96GV9', MSG
        assert res['RefSeq entry 1']['RefSeq ID'] == ['NP_149988.1'], MSG
        assert res['RefSeq entry 1']['sequence'] == seq, MSG
        assert res['RefSeq entry 1']['sequence matches'] == ['UniProt entry 1'], MSG
        assert res['RefSeq entry 1']['sequence length'] == 206, MSG

    def test_compare_entries_uniprotid(self) -> None:
        """Tests the comparison of database entries using a given example from a successful lookup by UniProt ID."""
        res = compare_entries(refseq_id=None, uniprot_id='Q9NY95', symbol=None)
        # data structure
        assert isinstance(res, pd.DataFrame)
        assert set(res.columns) == {'UniProt entry 1', 'RefSeq entry 1', 'RefSeq entry 2'}
        assert set(res.index) == {'symbol', 'UniProt ID', 'RefSeq ID', 'sequence', 'sequence matches',
                                  'sequence length'}
        # sequences
        seq_up1_rs2 = ('MEPPYSLTAHYDEFQEVKYVSRCGAGGARGASLPPGFPLGAARSATGARSGLPRWNRREVCLLSGLVFAAGLCAILAAMLALKYLGPVAAGGGAC'
                       + 'PEGCPERKAFARAARFLAANLDASIDPCQDFYSFACGGWLRRHAIPDDKLTYGTIAAIGEQNEERLRRLLARPGGGPGGAAQRKVRAFFRSCL'
                       + 'DMREIERLGPRPMLEVIEDCGGWDLGGAEERPGVAARWDLNRLLYKAQGVYSAAALFSLTVSLDDRNSSRYVIRIDQDGLTLPERTLYLAQDE'
                       + 'DSEKILAAYRVFMERVLSLLGADAVEQKAQEILQVEQQLANITVSEHDDLRRDVSSMYNKVTLGQLQKITPHLRWKWLLDQIFQEDFSEEEEV'
                       + 'VLLATDYMQQVSQLIRSTPHRVLHNYLVWRVVVVLSEHLSPPFREALHELAQEMEGSDKPQELARVCLGQANRHFGMALGALFVHEHFSAASK'
                       + 'AKVQQLVEDIKYILGQRLEELDWMDAETRAAARAKLQYMMVMVGYPDFLLKPDAVDKEYEFEVHEKTYFKNILNSIRFSIQLSVKKIRQEVDK'
                       + 'STWLLPPQALNAYYLPNKNQMVFPAGILQPTLYDPDFPQSLNYGGIGTIIGHELTHGYDDWGGQYDRSGNLLHWWTEASYSRFLRKAECIVRL'
                       + 'YDNFTVYNQRVNGKHTLGENIADMGGLKLAYHAYQKWVREHGPEHPLPRLKYTHDQLFFIAFAQNWCIKRRSQSIYLQVLTDKHAPEHYRVLG'
                       + 'SVSQFEEFGRAFHCPKDSPMNPAHKCSVW')
        seq_rs1 = ('MEPPYSLTAHYDEFQEVKYVSRCGAGGARGASLPPGFPLGAARSATGARSGLPRWNRREVCLLSGLVFAAGLCAILAAMLALKYLGPVAAGGGACPEGC'
                   + 'PERKAFARAARFLAANLDASIDPCQDFYSFACGGWLRRHAIPDDKLTYGTIAAIGEQNEERLRRLLARPGGGPGGAAQRKVRAFFRSCLDMREIERL'
                   + 'GPRPMLEVIEDCGGWDLGGAEERPGVAARWDLNRLLYKAQGVYSAAALFSLTVSLDDRNSSRYVIRIDQDGLTLPERTLYLAQDEDSEKILAAYRVF'
                   + 'MERVLSLLGADAVEQKAQEILQVEQQLANITVSEHDDLRRDVSSMYNKVTLGQLQKITPHLRWKWLLDQIFQEDFSEEEEVVLLATDYMQQVSQLIR'
                   + 'STPHRVLHNYLVWRVVVVLSEHLSPPFREALHELAQEMEGSDKPQELARVCLGQANRHFGMALGALFVHEHFSAASKAKVQQLVEDIKYILGQRLEE'
                   + 'LDWMDAETRAAARAKLQYMMVMVGYPDFLLKPDAVDKEYEFEVHEKTYFKNILNSIRFSIQLSVKKIRQEVDKWLLPPQALNAYYLPNKNQMVFPAG'
                   + 'ILQPTLYDPDFPQSLNYGGIGTIIGHELTHGYDDWGGQYDRSGNLLHWWTEASYSRFLRKAECIVRLYDNFTVYNQRVNGKHTLGENIADMGGLKLA'
                   + 'YHAYQKWVREHGPEHPLPRLKYTHDQLFFIAFAQNWCIKRRSQSIYLQVLTDKHAPEHYRVLGSVSQFEEFGRAFHCPKDSPMNPAHKCSVW')
        # values UniProt
        assert set(res['UniProt entry 1']['symbol']) == {'ECEL1', 'XCE'}, MSG
        assert res['UniProt entry 1']['UniProt ID'] == 'Q9NY95', MSG
        assert set(res['UniProt entry 1']['RefSeq ID']) == {'NP_001277716.1', 'NP_004817.2'}, MSG
        assert res['UniProt entry 1']['sequence'] == seq_up1_rs2, MSG
        assert res['UniProt entry 1']['sequence matches'] == ['RefSeq entry 2'], MSG
        assert res['UniProt entry 1']['sequence length'] == 775, MSG
        # values RefSeq 1
        assert res['RefSeq entry 1']['symbol'] == ['ECEL1'], MSG
        assert res['RefSeq entry 1']['UniProt ID'] == 'O95672', MSG
        assert res['RefSeq entry 1']['RefSeq ID'] == ['NP_001277716.1'], MSG
        assert res['RefSeq entry 1']['sequence'] == seq_rs1, MSG
        assert not res['RefSeq entry 1']['sequence matches'], MSG
        assert res['RefSeq entry 1']['sequence length'] == 773, MSG
        # values RefSeq 2
        assert res['RefSeq entry 2']['symbol'] == ['ECEL1'], MSG
        assert res['RefSeq entry 2']['UniProt ID'] == 'O95672', MSG
        assert res['RefSeq entry 2']['RefSeq ID'] == ['NP_004817.2'], MSG
        assert res['RefSeq entry 2']['sequence'] == seq_up1_rs2, MSG
        assert res['RefSeq entry 2']['sequence matches'] == ['UniProt entry 1'], MSG
        assert res['RefSeq entry 2']['sequence length'] == 775, MSG

    def test_compare_entries_symbol(self) -> None:
        """Tests the comparison of database entries using a given example from a query by symbol."""
        res = compare_entries(refseq_id=None, uniprot_id=None, symbol='OXT')
        # data structure
        assert isinstance(res, pd.DataFrame)
        assert set(res.columns) == {'UniProt entry 1', 'No RefSeq entry'}, MSG
        assert set(res.index) == {'symbol', 'UniProt ID', 'RefSeq ID', 'sequence', 'sequence matches',
                                  'sequence length'}
        # values UniProt
        assert set(res['UniProt entry 1']['symbol']) == {'OXT', 'OT'}, MSG
        assert res['UniProt entry 1']['UniProt ID'] == 'Q3MIG0', MSG
        assert set(res['UniProt entry 1']['RefSeq ID']) == {'NP_000906.1', 'XP_011527540.1'}, MSG
        assert res['UniProt entry 1']['sequence'] == ('MAGPSLACCLLGLLALTSACYIQNCPLGGKRAAPDLDVRKCLPCGPGGKGRCFGPNICCA'
               + 'EELGCFVGTAEALRCQEENYLPSPCQSGQKACGSGGRCAVLGLCCSPDGCHADPACDAEATFSQR'), MSG
        assert not res['UniProt entry 1']['sequence matches'], MSG
        assert res['UniProt entry 1']['sequence length'] == 125, MSG
        # values RefSeq
        assert not any(res['No RefSeq entry'].values), MSG

    def test_summary_statistics(self) -> None:
        """Determines the correctness of calculated summary statistics in a range of +/- 5% of a comparison value."""
        res = summary_statistics()
        # data structure
        assert isinstance(res, pd.DataFrame)
        assert set(res.columns) == {'Number of matches', 'Matching UniProt entries [%]', 'Matching RefSeq entries [%]'}
        assert set(res.index) == {'Symbol', 'RefSeq ID', 'UniProt ID', 'Sequence', 'Sequence length'}
        # values - ballpark figures
        matches = res['Number of matches']
        assert is_approx(matches['Symbol'], 16281), MSG
        assert is_approx(matches['RefSeq ID'], 39343), MSG
        assert is_approx(matches['UniProt ID'], 1802), MSG
        assert is_approx(matches['Sequence'], 22095), MSG
        assert is_approx(matches['Sequence length'], 22525), MSG
        uniprot = res['Matching UniProt entries [%]']
        assert is_approx_percentage(uniprot['Symbol'], 32.26), MSG
        assert is_approx_percentage(uniprot['RefSeq ID'], 89.32), MSG
        assert is_approx_percentage(uniprot['UniProt ID'], 7.07), MSG
        assert is_approx_percentage(uniprot['Sequence'], 85.48), MSG
        assert is_approx_percentage(uniprot['Sequence length'], 86.63), MSG
        refseq = res['Matching RefSeq entries [%]']
        assert is_approx_percentage(refseq['Symbol'], 26.67), MSG
        assert is_approx_percentage(refseq['RefSeq ID'], 64.45), MSG
        assert is_approx_percentage(refseq['UniProt ID'], 2.95), MSG
        assert is_approx_percentage(refseq['Sequence'], 36.19), MSG
        assert is_approx_percentage(refseq['Sequence length'], 36.9), MSG

    def test_update_stats(self) -> None:
        """Tests summary_statistics' helper function with two consecutive calls on a dummy consensus collection."""
        # empty collection
        consensus = {key: {"matches": 0, "UniProt entry": set()}
                     for key in ["Symbol", "RefSeq ID", "UniProt ID", "Sequence", "Sequence length"]}
        # query for actin beta from UniProt
        refseq_id = 'NP_001092.1'
        refseq_entry = get_refseq_entry(refseq_id)
        uniprot_id = 'Q96HG5'
        uniprot_entry = get_uniprot_entry(uniprot_id)
        res = update_stats(refseq_entry, uniprot_entry, refseq_id, uniprot_id, "UniProt", consensus)
        assert res['Symbol']['matches'] == 1
        assert res['Symbol']['UniProt entry'] == {uniprot_id}
        assert res['RefSeq ID']['matches'] == 0
        assert not res['RefSeq ID']['UniProt entry']
        assert res['UniProt ID']['matches'] == 0
        assert not res['UniProt ID']['UniProt entry']
        assert res['Sequence']['matches'] == 1
        assert res['Sequence']['UniProt entry'] == {uniprot_id}
        assert res['Sequence length']['matches'] == 1
        assert res['Sequence length']['UniProt entry'] == {uniprot_id}
        # query from RefSeq
        refseq_id = 'NP_001303897.1'
        refseq_entry = get_refseq_entry(refseq_id)
        uniprot_id = 'Q96GV9'
        uniprot_entry = get_uniprot_entry(uniprot_id)
        res = update_stats(refseq_entry, uniprot_entry, refseq_id, uniprot_id, "RefSeq", consensus)
        assert res['Symbol']['matches'] == 1
        assert res['Symbol']['UniProt entry'] == {'Q96HG5'}
        assert res['RefSeq ID']['matches'] == 0
        assert not res['RefSeq ID']['UniProt entry']
        assert res['UniProt ID']['matches'] == 0
        assert not res['UniProt ID']['UniProt entry']
        assert res['Sequence']['matches'] == 2
        assert res['Sequence']['UniProt entry'] == {'Q96GV9', 'Q96HG5'}
        assert res['Sequence length']['matches'] == 2
        assert res['Sequence length']['UniProt entry'] == {'Q96GV9', 'Q96HG5'}

    def test_finalize_stats(self) -> None:
        """Checks the correct conversion of collected consensus between databases as a dictionary to a pd.DataFrame."""
        consensus = {'Symbol': {'matches': 12, 'UniProt entry': set('ABCDEFGHIJK')},
                     'RefSeq ID': {'matches': 6, 'UniProt entry': set('ABCDEF')},
                     'UniProt ID': {'matches': 6, 'UniProt entry': set('ABCDE')},
                     'Sequence': {'matches': 2, 'UniProt entry': set('AB')},
                     'Sequence length': {'matches': 3, 'UniProt entry': set('AB')}}
        res = finalize_stats(consensus, 21, 32)
        # data structure
        assert isinstance(res, pd.DataFrame)
        assert set(res.columns) == {'Number of matches', 'Matching UniProt entries [%]', 'Matching RefSeq entries [%]'}
        assert set(res.index) == {'Symbol', 'RefSeq ID', 'UniProt ID', 'Sequence', 'Sequence length'}
        # values
        matches = res['Number of matches']
        assert matches['Symbol'] == 12
        assert matches['RefSeq ID'] == 6
        assert matches['UniProt ID'] == 6
        assert matches['Sequence'] == 2
        assert matches['Sequence length'] == 3
        uniprot = res['Matching UniProt entries [%]']
        assert uniprot['Symbol'] == '52.38%'
        assert uniprot['RefSeq ID'] == '28.57%'
        assert uniprot['UniProt ID'] == '23.81%'
        assert uniprot['Sequence'] == '9.52%'
        assert uniprot['Sequence length'] == '9.52%'
        refseq = res['Matching RefSeq entries [%]']
        assert refseq['Symbol'] == '37.50%'
        assert refseq['RefSeq ID'] == '18.75%'
        assert refseq['UniProt ID'] == '18.75%'
        assert refseq['Sequence'] == '6.25%'
        assert refseq['Sequence length'] == '9.38%'
