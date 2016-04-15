import mmt_defaults
import json_config
import os
from data_classes import Datafile
from tasks_v2 import Task
from external_tools import TOOLS_DICT


busco_defaults = {'arthropoda': False, 'metazoa': False,
                  'vertebrata': False, 'eukaryota': False,
                  'fungi': False, 'bacteria': False,
                  'plantae': False}



def gen_db_supervisor(force=False, sprot=False, uniref90=False nr=False, busco_args=busco_defaults, blast_plus=False, cpu=float('inf'), dep):
    db_configs = get_db_config()
    tasks = []
    if(sprot):
        sprot_dep = []
        if(force or not os.path.exists(db_configs['uniprot_sprot_fasta'])):
            sprot_download = download_task(mmt_defaults.URL_UNIPROT_SPROT, db_configs['uniprot_sprot_fasta'], '', [])
            tasks.append(sprot_download)
            sprot_dep.append(sprot_download)
        sprot_install = 



if(__name__ == '__main__'):
    pass
