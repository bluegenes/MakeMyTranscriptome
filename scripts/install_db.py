import os
from data_classes import get_dbs
from tasks_v2 import Supervisor
import functions_databases as fdb
import argparse
import mmt_defaults as statics

busco_defaults = {'arthropoda': False, 'metazoa': False,
                  'vertebrata': False, 'eukaryota': False,
                  'fungi': False, 'bacteria': False,
                  'plantae': False}


def download_task_wrapper(db, tasks):
    return fdb.download_task(db.url, db.download_location, db.type, tasks)


def gen_dmnd_blast_tasks(db, force, blast_plus):
    tasks = []
    depends = []
    sprot_download = download_task_wrapper(db, [])
    tasks.append(sprot_download)
    depends.append(sprot_download)
    install_dmnd = fdb.build_diamond_task(db.download_location, db.call_path, depends)
    tasks.append(install_dmnd)
    if(blast_plus):
        install_blast = fdb.build_blast_task(db.download_location, db.call_path, 'prot', depends)
        tasks.append(install_blast)
    return Supervisor(tasks)


def check_db_dir():
    def make_dir(_path):
        if(not os.path.exists(_path)):
            os.mkdir(_path)

    dirs = [statics.PATH_DATABASE_LOGS,
            statics.PATH_UNIPROT_SPROT_DIR,
            statics.PATH_UNIREF90_DIR,
            statics.PATH_NR_DIR,
            statics.PATH_PFAM_DIR,
            statics.PATH_BUSCO_REFERENCE_DIR]
    for d in dirs:
        make_dir(d)


def gen_db_supervisor(force=False, sprot=False, uniref90=False, nr=False, busco_args=busco_defaults, blast_plus=False, idmapping=False, cpu=float('inf'), dep=[]):
    check_db_dir()
    dbs = get_dbs(defaults=force)
    tasks = []
    if(sprot):
        tasks.append(gen_dmnd_blast_tasks(dbs['uniprot_sprot'], force, blast_plus))
    if(uniref90):
        tasks.append(gen_dmnd_blast_tasks(dbs['uniref90'], force, blast_plus))
    if(nr):
        tasks.append(gen_dmnd_blast_tasks(dbs['nr'], force, blast_plus))
    for busco_db in busco_args:
        if(busco_args[busco_db]):
            tasks.append(download_task_wrapper(dbs['busco_'+busco_db], []))
    if(idmapping):
        tasks.append(download_task_wrapper(dbs['id_mapping'], []))
        tasks.append(download_task_wrapper(dbs['idmapping_selected'], []))
    for db_string in dbs:
        if(db_string in ['id_mapping', 'idmapping_selected','uniprot_sprot', 'uniref90', 'nr'] or db_string.startswith('busco_')):
            pass
        else:
            tasks.append(download_task_wrapper(dbs[db_string], []))
    tasks = [t for t in tasks if(t is not None)]
    return Supervisor(tasks, cpu=cpu)


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument('--hard', action='store_true', default=False)
    parser.add_argument('--swiss_prot', action='store_true', default=False)
    parser.add_argument('--uniref90', action='store_true', default=False)
    parser.add_argument('--nr', action='store_true', default=False)
    parser.add_argument('--buscos', help='a comma seperated list of busco files that need to be downloaded')
    parser.add_argument('--cpu', type=int, default=4)
    parser.add_argument('--buildBlastPlus', action='store_true', default=False)
    parser.add_argument('--idmapping', action='store_true', default=False)
    args = parser.parse_args()
    busco_flags = {}
    if(args.buscos is not None):
        args.buscos = args.buscos.split(',')
        for b in args.buscos:
            busco_flags[b] = True
    sup = gen_db_supervisor(args.hard, args.swiss_prot, args.uniref90, args.nr, busco_flags, args.buildBlastPlus, args.idmapping, args.cpu)
    sup.run()

