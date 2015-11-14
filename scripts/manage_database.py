from task_functions_v2 import (
    PATH_DATABASES, PATH_UNIREF90, PATH_SWISS_PROT,
    PATH_NR, PATH_BUSCO_METAZOA, pfam_build_task,
    build_diaimond_task, build_blast_task)
from time import strftime
import gzip
import tarfile
import os
import json
import sys
from tasks_v2 import Supervisor
if(sys.version[0] == '3'):
    from urllib.request import urlretrieve, ContentTooShortError
else:
    from urllib import urlretrieve, ContentTooShortError

swissprot_folder = os.path.join(PATH_DATABASES, 'uniprot_sprot')
uniref90_folder = os.path.join(PATH_DATABASES, 'uniref90')
nr_folder = os.path.join(PATH_DATABASES, 'nr')

'''
makeblastdb *3
make daiomond *3
pfam
'''

url_sprot = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
sprot_target = PATH_SWISS_PROT

url_uniref90 = 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz'
uniref90_target = PATH_UNIREF90

url_nr = 'ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz'
nr_target = PATH_NR

url_busco_metazoa = 'http://busco.ezlab.org/files/metazoa_buscos.tar.gz'
busco_metazoa_target = PATH_BUSCO_METAZOA

url_go_pathway = 'http://rest.genome.jp/link/go/pathway'
go_pathway_target = '{0!s}/go_pathway.txt'.format(PATH_DATABASES)

url_swiss_enzyme = 'http://rest.genome.jp/link/enzyme/swissprot'
swiss_enzyme_target = '{0!s}/swiss_enzyme.list'.format(PATH_DATABASES)

url_pfam_enzyme = 'http://rest.genome.jp/link/enzyme/pfam'
pfam_enzyme_target = '{0!s}/pfam_enzyme.list'.format(PATH_DATABASES)

url_id_mapping = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz'
id_mapping_target = '{0!s}/idmapping.dat'.format(PATH_DATABASES)

url_sprot_map = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
sprot_map_target = '{0!s}/uniprot_sprot.dat'.format(PATH_DATABASES)

url_kog_functional = 'http://eggnogdb.embl.de/download/eggnog_4.1/data/NOG/NOG.annotations.tsv.gz'
kog_functional_target = '{0!s}/allKOG_functional_info.txt'.format(PATH_DATABASES)

url_slim_generic = 'http://www.geneontology.org/ontology/subsets/goslim_generic.obo'
slim_generic_target = '{0!s}/goslim_generic.obo'.format(PATH_DATABASES)

url_pfam_db = 'https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/Pfam-A.hmm.gz'
pfam_db_target = '{0!s}/Pfam-A.hmm'.format(PATH_DATABASES)

url_idmapping_selected = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz'
idmapping_selected_target = '{0!s}/idmapping_selected.tab'.format(PATH_DATABASES)

idmapping_keys = {'BioCyc': '{0!s}/idmapping.biocyc'.format(PATH_DATABASES),
                  'eggNOG': '{0!s}/idmapping.eggNOG'.format(PATH_DATABASES),
                  'KO': '{0!s}/idmapping.KO'.format(PATH_DATABASES),
                  'OrthoDB': '{0!s}/idmapping.orthodb'.format(PATH_DATABASES)}


def run_tasks(tasks):
    for t in tasks:
        t.stdout = t.name+'.stdout'
        t.stderr = t.name+'.stderr'

    s = Supervisor(tasks=tasks, force_run=True)
    s.run()
    for t in tasks:
        os.remove(t.stdout)
        os.remove(t.stderr)


def safe_retrieve(source, target):
    print('getting '+source)
    urlretrieve(source, target+'.temp')
    os.rename(target+'.temp', target)


def url_unzip(source, target):
    print('getting '+source)
    urlretrieve(source, target+'.gz')
    f = gzip.open(target+'.gz', 'rb')
    g = open(target, 'wb')
    for line in f:
        g.write(line)
    f.close()
    g.close()
    os.remove(target+'.gz')


def tar_retrieve(source, target):
    print('getting '+source)
    urlretrieve(source, target+'.tar.gz')
    tfile = tarfile.open(source, 'r:gz')
    tfile.extractall(target)
    os.remove('.tar.gz')


def get(log_table, flag, source, target):
    try:
        if(flag == 'gz'):
            url_unzip(source, target)
        elif(flag == ''):
            safe_retrieve(source, target)
        elif(flag == 'tar'):
            tar_retrieve(source, target)
        else:
            print('Your shits fucked.')
        basename = os.path.basename(target)
        log_table[basename] = strftime('%b-%d-%Y')
    except ContentTooShortError:
        print('failed to install {0!s}'.format(basename(arget)))


def read_log():
    log = open(PATH_DATABASES+'/database_log')
    log_table = json.load(log)
    log.close()
    return log_table


def write_log(log_table):
    log = open(PATH_DATABASES+'/database_log', 'w')
    json.dump(log_table, log, sort_keys=True, indent=4)
    log.close()


def download_databases(log_table):
    partial_get = lambda a, b, c: get(log_table, a, b, c)
    partial_get('', url_go_pathway, go_pathway_target)
    partial_get('', url_swiss_enzyme, swiss_enzyme_target)
    partial_get('', url_pfam_enzyme, pfam_enzyme_target)
    partial_get('', url_slim_generic, slim_generic_target)
    partial_get('gz', url_sprot, sprot_target)
    partial_get('gz', url_uniref90, uniref90_target)
    partial_get('gz', url_nr, nr_target)
    partial_get('gz', url_id_mapping, id_mapping_target)
    partial_get('gz', url_kog_functional, kog_functional_target)
    partial_get('gz', url_pfam_db, pfam_db_target)
    partial_get('tar', url_busco_metazoa, busco_metazoa_target)
    return log_table


def subset_dat(dat_file, key_file_dict, log_table):
    date = log_table[basename(dat_file)]
    for key in key_file_dict:
        key = basename(key_file_dict[key])
        log_table[key] = date
    dat = open(dat_file)
    key_file = {key: open(key_file_dict[key], 'w') for key in key_file_dict}
    for line in dat:
        values = line.split('\t')
        if(values[1] in key_file):
            key_file[values[1]].write(line)
    for key in key_file:
        key_file[key].close()
    return log_table


def check_database_dir():
    if(not os.path.isdir(PATH_DATABASES)):
        os.mkdir(PATH_DATABASES)
        write_log({})
    if(not os.path.isdir(swissprot_folder)):
        os.mkdir(swissprot_folder)
    if(not os.path.isdir(uniref90_folder)):
        os.mkdir(uniref90_folder)
    if(not os.path.isdir(nr_folder)):
        os.mkdir(nr_folder)


def main():
    check_database_dir()
    log_table = read_log()
    log_table = download_databases(log_table)
    log_table = subset_dat(id_mapping_target, idmapping_keys, log_table)
    nr_task = build_blast_task(nr_target, nr_target, 'prot', [])
    nr_diamond = build_diaimond_task(nr_target, nr_target, [])
    swissprot_task = build_blast_task(sprot_target, sprot_target, 'prot', [])
    swissprot_diamond = build_diaimond_task(sprot_target, sprot_target, [])
    uniref90_task = build_blast_task(uniref90_target, uniref90_target, 'prot', [])
    uniref90_diamond = build_diaimond_task(uniref90_target, uniref90_target, [])
    pfam_task = pfam_build_task(pfam_db_target, [])
    run_tasks([nr_task, nr_diamond, swissprot_task, swissprot_diamond,
               uniref90_task, uniref90_diamond, pfam_task])
    write_log(log_table)


if(__name__ == '__main__'):
    main()
