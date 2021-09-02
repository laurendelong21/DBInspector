from pathlib import Path
import os
import os.path as osp
import logging

# set up directories
HOME = str(Path.home())
CACHE = osp.join(HOME, '.dbinspector')
LOGS = osp.join(CACHE, 'logs')
# parsed
UNIPROT = osp.join(CACHE, 'uniprot')
REFSEQ = osp.join(CACHE, 'refseq')
# downloads
DATA = osp.join(CACHE, 'data')
REFSEQ_FASTA = osp.join(DATA, 'refseq_fasta')

for folder in [CACHE, LOGS, UNIPROT, REFSEQ, DATA, REFSEQ_FASTA]:
    os.makedirs(folder, exist_ok=True)

# logging
logging.basicConfig(filename=osp.join(LOGS, 'dbinspection.log'),
                    filemode='a',
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%d/%m/%Y %I:%M:%S %p')
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
