import os.path as osp
from dbinspector.startup import REFSEQ, UNIPROT

from dbinspector.map import read_refseq_data, read_uniprot_data, get_refseq_entry, get_uniprot_entry, \
    retrieve_by_refseq_id, retrieve_by_uniprot_id, retrieve_by_symbol, find_entries
from dbinspector.parse import parse_all

MSG = "This test might be compromised by updated metadata from the database."
NONSENSE = 'Blubb~Blubb~Blubb~O0o.'


def metadata_keys_complete(result_keys) -> bool:
    """
    Helper function for tests that assert, that all metdata keys are present in a dictionary that has been
    returned as the result of a query.

    :param result_keys: keys of a database-specific dictionary that has been returned by a query, corresponds to
                        metadata categories
    :type result_keys: dict_keys
    """
    return set(result_keys) == {'symbol', 'UniProt ID', 'sequence', 'RefSeq ID'}


class TestMap:
    """Test class for the map module"""

    def test_preparation(self):
        """Makes sure that the parsed data is available which is the underlying assumptions of functions in map."""
        parse_all()
        refseq_path = osp.join(REFSEQ, 'refseq.json')
        uniprot_path = osp.join(UNIPROT, 'uniprot.json')
        assert osp.isfile(refseq_path)
        assert osp.isfile(uniprot_path)

    def test_read_parsed_data(self):
        """Checks that parsed data is loaded/read in correctly from the cached json files to a dictionary."""
        refseq = read_refseq_data()
        uniprot = read_uniprot_data()
        # data structure
        assert isinstance(refseq, dict)
        assert all([isinstance(key, str) for key in refseq.keys()])
        assert all([isinstance(val, dict) for val in refseq.values()])
        assert isinstance(uniprot, dict)
        assert all([isinstance(key, str) for key in uniprot.keys()])
        assert all([isinstance(val, dict) for val in uniprot.values()])
        # all refseq keys (= acc. id) is NP_...
        assert all([k.startswith('NP_') for k in refseq.keys()])
        # number of entries
        assert 0.95*61047 < len(refseq) < 1.05*61047, MSG
        assert 0.95*20312 < len(uniprot) < 1.05*20312, MSG

    def test_get_refseq_entry(self):
        """Tests the lookup of a RefSeq entry by accession ID by comparison to an example."""
        # nonsense
        assert not get_refseq_entry(NONSENSE)
        # successful
        res: dict = get_refseq_entry('NP_149988.1')
        assert metadata_keys_complete(res.keys())
        assert not res['symbol'], MSG
        assert res['UniProt ID'] == 'Q96GV9', MSG
        assert res['sequence'] == ('MEVDINGESRSTLTTLPFPGAEANSPGKAEAEKPRCSSTPCSPMRRTVSGYQILHMDSNYLVGFTTGEELLKLAQKCTGGEES'
               + 'KAEAMPSLRSKQLDAGLARSSRLYKTRSRYYQPYEIPAVNGRRRRRMPSSGDKCTKSLPYEPYKALHGPLPLCLLKGKRAHSKSLDYLNLDKMIKEPADTE'
               + 'VLQYQLQHLTLRGDRVFARNNT'), MSG
        assert res['RefSeq ID'] == ['NP_149988.1'], MSG

    def test_get_uniprot_entry(self):
        """Tests the lookup of a UniProt entry by accession ID by comparison to an example."""
        # nonsense
        assert not get_uniprot_entry(NONSENSE)
        # successful
        res: dict = get_uniprot_entry('Q9NSV2')
        assert metadata_keys_complete(res.keys())
        assert set(res['symbol']) == {'GOLGA2P5', 'GOLGA2B', 'GOLGA2L1'}, MSG
        assert res['UniProt ID'] == 'Q9NSV2', MSG
        assert res['sequence'] == ('MDSEEEEEVPQPMPSIPEDLESQKAMVAFFNSAVASAEEEQARLCGQLKECTASAWLICWPRPRRNLRQQPQPQELGVIPCVG'
                                   + 'RPTRPCRGPWRSCGRVHRTVPEPEGSAEGGGVHQQAGPGQGRGEGEAAGAGVACGRLQQVA'), MSG
        assert not res['RefSeq ID'], MSG

    def test_retrieve_by_refseq_id(self):
        """Tests retrieval of appropriate RefSeq and UniProt entries upon query by RefSeq accession ID, checks
        correctness by comparison to an example."""
        # ==== nonsense ====
        assert not any(retrieve_by_refseq_id(NONSENSE).values())
        # ==== successful ==== 1 RefSeq & 1 UniProt entry
        res: dict = retrieve_by_refseq_id('NP_149988.1')
        assert res['RefSeq'] and res['UniProt']
        # RefSeq
        assert isinstance(res['RefSeq'], list)
        assert len(res['RefSeq']) == 1, MSG
        assert isinstance(res['RefSeq'][0], dict), MSG
        assert metadata_keys_complete(res['RefSeq'][0].keys())
        assert not res['RefSeq'][0]['symbol'], MSG
        assert res['RefSeq'][0]['UniProt ID'] == 'Q96GV9', MSG
        assert res['RefSeq'][0]['sequence'] == ('MEVDINGESRSTLTTLPFPGAEANSPGKAEAEKPRCSSTPCSPMRRTVSGYQILHMDSNYLVGFTTGEEL'
               + 'LKLAQKCTGGEESKAEAMPSLRSKQLDAGLARSSRLYKTRSRYYQPYEIPAVNGRRRRRMPSSGDKCTKSLPYEPYKALHGPLPLCLLKGKRAHSKSLDYL'
               + 'NLDKMIKEPADTEVLQYQLQHLTLRGDRVFARNNT'), MSG
        assert isinstance(res['RefSeq'][0]['RefSeq ID'], list), MSG
        assert len(res['RefSeq'][0]['RefSeq ID']) == 1, MSG
        assert res['RefSeq'][0]['RefSeq ID'][0] == 'NP_149988.1', MSG
        # UniProt
        assert isinstance(res['UniProt'], list)
        assert len(res['UniProt']) == 1, MSG
        assert isinstance(res['UniProt'][0], dict), MSG
        assert metadata_keys_complete(res['UniProt'][0].keys())
        assert set(res['UniProt'][0]['symbol']) == {'MACIR', 'C5orf30'}, MSG
        assert res['UniProt'][0]['UniProt ID'] == 'Q96GV9', MSG
        assert res['UniProt'][0]['sequence'] == ('MEVDINGESRSTLTTLPFPGAEANSPGKAEAEKPRCSSTPCSPMRRTVSGYQILHMDSNYLVGFTTGEE'
               + 'LLKLAQKCTGGEESKAEAMPSLRSKQLDAGLARSSRLYKTRSRYYQPYEIPAVNGRRRRRMPSSGDKCTKSLPYEPYKALHGPLPLCLLKGKRAHSKSLDY'
               + 'LNLDKMIKEPADTEVLQYQLQHLTLRGDRVFARNNT'), MSG
        assert isinstance(res['UniProt'][0]['RefSeq ID'], list), MSG
        assert set(res['UniProt'][0]['RefSeq ID']) == {'NP_001303897.1', 'NP_001303898.1', 'NP_149988.1'}, MSG

    def test_retrieve_by_uniprot_id(self):
        """Tests retrieval of appropriate RefSeq and UniProt entries upon query by RefSeq accession ID, checks
        correctness by comparison to an example without corresponding RefSeq entry."""
        # ==== nonsense ====
        assert not any(retrieve_by_uniprot_id(NONSENSE).values())
        # ==== successful ==== NO RefSeq & 1 UniProt entry
        res: dict = retrieve_by_uniprot_id('Q9NSV2')
        assert res['UniProt'] and not res['RefSeq']
        assert isinstance(res['UniProt'], list)
        assert len(res['UniProt']) == 1, MSG
        assert isinstance(res['UniProt'][0], dict), MSG
        assert metadata_keys_complete(res['UniProt'][0].keys())
        assert set(res['UniProt'][0]['symbol']) == {'GOLGA2P5', 'GOLGA2B', 'GOLGA2L1'}, MSG
        assert res['UniProt'][0]['UniProt ID'] == 'Q9NSV2', MSG
        assert res['UniProt'][0]['sequence'] == ('MDSEEEEEVPQPMPSIPEDLESQKAMVAFFNSAVASAEEEQARLCGQLKECTASAWLICWPRPRRNLRQ'
               + 'QPQPQELGVIPCVGRPTRPCRGPWRSCGRVHRTVPEPEGSAEGGGVHQQAGPGQGRGEGEAAGAGVACGRLQQVA'), MSG
        assert not res['UniProt'][0]['RefSeq ID'], MSG

    def test_retrieve_by_symbol(self):
        """Tests that retrieval by symbol returns entries from UniProt and RefSeq in the appropriate data structure."""
        # ==== nonsense ====
        assert not any(retrieve_by_symbol(NONSENSE).values())
        # ==== successful ==== actin beta ==== 1 RefSeq & 1 UniProt entry
        res: dict = retrieve_by_symbol('ACTB')
        assert res['RefSeq'] and res['UniProt']
        actinb_seq = ('MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRV'
                      + 'APEEHPVLLTEAPLNPKANREKMTQIMFETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDGVTHTVPIYEGYALPHAILRLDLAGRDLTDYLM'
                      + 'KILTERGYSFTTTAEREIVRDIKEKLCYVALDFEQEMATAASSSSLEKSYELPDGQVITIGNERFRCPEALFQPSFLGMESCGIHETTFNSIMK'
                      + 'CDVDIRKDLYANTVLSGGTTMYPGIADRMQKEITALAPSTMKIKIIAPPERKYSVWIGGSILASLSTFQQMWISKQEYDESGPSIVHRKCF')
        # RefSeq
        assert isinstance(res['RefSeq'], list)
        assert len(res['RefSeq']) == 1, MSG
        assert isinstance(res['RefSeq'][0], dict), MSG
        assert metadata_keys_complete(res['RefSeq'][0].keys())
        assert res['RefSeq'][0]['symbol'] == ['ACTB'], MSG
        assert res['RefSeq'][0]['UniProt ID'] == 'Q1KLZ0', MSG
        assert res['RefSeq'][0]['sequence'] == actinb_seq, MSG
        assert res['RefSeq'][0]['RefSeq ID'] == 'NP_001092.1', MSG
        # UniProt
        assert isinstance(res['UniProt'], list)
        assert len(res['UniProt']) == 1, MSG
        assert isinstance(res['UniProt'][0], dict), MSG
        assert metadata_keys_complete(res['UniProt'][0].keys())
        assert res['UniProt'][0]['symbol'] == ['ACTB'], MSG
        assert res['UniProt'][0]['UniProt ID'] == 'Q96HG5', MSG
        assert res['UniProt'][0]['sequence'] == actinb_seq, MSG
        assert res['UniProt'][0]['RefSeq ID'] == ['NP_001092.1'], MSG
        # ==== successful ==== Phenylethanolamine N-methyltransferase ==== 1 RefSeq & 2 UniProt entries
        res: dict = retrieve_by_symbol('PNMT')
        assert res['RefSeq'] and res['UniProt']
        # RefSeq
        assert isinstance(res['RefSeq'], list)
        assert len(res['RefSeq']) == 1, MSG
        assert isinstance(res['RefSeq'][0], dict), MSG
        assert metadata_keys_complete(res['RefSeq'][0].keys())
        assert res['RefSeq'][0]['symbol'] == ['PNMT'], MSG
        assert res['RefSeq'][0]['UniProt ID'] == 'P11086', MSG
        assert res['RefSeq'][0]['sequence'] == ('MSGADRSPNAGAAPDSAPGQAAVASAYQRFEPRAYLRNNYAPPRGDLCNPNGVGPWKLRCLAQTFATGEV'
               + 'SGRTLIDIGSGPTVYQLLSACSHFEDITMTDFLEVNRQELGRWLQEEPGAFNWSMYSQHACLIEGKGECWQDKERQLRARVKRVLPIDVHQPQPLGAGSPA'
               + 'PLPADALVSAFCLEAVSPDLASFQRALDHITTLLRPGGHLLLIGALEESWYLAGEARLTVVPVSEEEVREALVRSGYKVRDLRTYIMPAHLQTGVDDVKGV'
               + 'FFAWAQKVGL'), MSG
        assert res['RefSeq'][0]['RefSeq ID'] == 'NP_002677.1', MSG
        # UniProt
        assert isinstance(res['UniProt'], list)
        assert len(res['UniProt']) == 2, MSG
        assert isinstance(res['UniProt'][0], dict), MSG
        assert metadata_keys_complete(res['UniProt'][0].keys())
        assert set(res['UniProt'][0]['symbol']) == {'PEMT', 'PEMPT', 'PNMT'}, MSG
        assert res['UniProt'][0]['UniProt ID'] == 'Q9Y6V9', MSG
        assert res['UniProt'][0]['sequence'] == ('MTRLLGYVDPLDPSFVAAVITITFNPLYWNVVARWEHKTRKLSRAFGSPYLACYSLSVTILLLNFLRSH'
               + 'CFTQAMLSQPRMESLDTPAAYSLGLALLGLGVVLVLSSFFALGFAGTFLGDYFGILKEARVTVFPFNILDNPMYWGSTANYLGWAIMHASPTGLLLTVLVA'
               + 'LTYIVALLYEEPFTAEIYRQKASGSHKRS'), MSG
        assert set(res['UniProt'][0]['RefSeq ID']) == {'NP_001254480.1', 'NP_001254481.1', 'NP_009100.2',
                                                       'NP_680477.1', 'NP_680478.1'}, MSG
        assert isinstance(res['UniProt'][1], dict), MSG
        assert metadata_keys_complete(res['UniProt'][1].keys())
        assert set(res['UniProt'][1]['symbol']) == {'PNMT', 'PENT'}, MSG
        assert res['UniProt'][1]['UniProt ID'] == 'P11086', MSG
        assert res['UniProt'][1]['sequence'] == ('MSGADRSPNAGAAPDSAPGQAAVASAYQRFEPRAYLRNNYAPPRGDLCNPNGVGPWKLRCLAQTFATGE'
               + 'VSGRTLIDIGSGPTVYQLLSACSHFEDITMTDFLEVNRQELGRWLQEEPGAFNWSMYSQHACLIEGKGECWQDKERQLRARVKRVLPIDVHQPQPLGAGSP'
               + 'APLPADALVSAFCLEAVSPDLASFQRALDHITTLLRPGGHLLLIGALEESWYLAGEARLTVVPVSEEEVREALVRSGYKVRDLRTYIMPAHLQTGVDDVKG'
               + 'VFFAWAQKVGL'), MSG
        assert res['UniProt'][1]['RefSeq ID'] == ['NP_002677.1'], MSG

    def test_find_entries(self):
        """Tests the appropriate functions' return data types."""
        assert find_entries(refseq_id='NP_149988.1')['RefSeq']
        assert find_entries(uniprot_id='Q9NSV2')['UniProt']
        # ==== nonsense ====
        assert not any(find_entries(refseq_id=NONSENSE).values())
        assert not any(find_entries(uniprot_id=NONSENSE).values())
        assert not any(find_entries(symbol=NONSENSE).values())
