from path_vars import PATH_DATABASES
from time import strftime
import gzip
from os import rename, remove
from os.path import basename
import json
import sys
if(sys.version[0] == '3'):
    from urllib.request import urlretrieve, ContentTooShortError
else:
    from urllib import urlretrieve, ContentTooShortError

url_sprot = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
sprot_target = '{0!s}/uniprot_sprot/uniprot_sprot.fasta'.format(PATH_DATABASES)

url_uniref90 = 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz'
uniref90_target = '{0!s}/uniref90/uniref90.fasta'.format(PATH_DATABASES)

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

url_idmapping_selected = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz'
idmapping_selected_target = '{0!s}/idmapping_selected.tab'.format(PATH_DATABASES)

idmapping_keys = {'BioCyc': '{0!s}/idmapping.biocyc'.format(PATH_DATABASES),
                  'eggNOG': '{0!s}/idmapping.eggNOG'.format(PATH_DATABASES),
                  'KO': '{0!s}/idmapping.KO'.format(PATH_DATABASES),
                  'OrthoDB': '{0!s}/idmapping.orthodb'.format(PATH_DATABASES)}


def safe_retrieve(source, target):
    print('getting '+source)
    urlretrieve(source, target+'.temp')
    rename(target+'.temp', target)


def url_unzip(source, target):
    print('getting '+source)
    urlretrieve(source, target+'.gz')
    f = gzip.open(target+'.gz', 'rb')
    g = open(target, 'wb')
    for line in f:
        g.write(line)
    f.close()
    g.close()
    remove(target+'.gz')


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
    date = strftime('%b-%d-%Y')

    try:
        if(not os.path.isdir(PATH_DATABASES+'/uniprot_sprot')):
            os.mkdirs(PATH_DATABASES+'/uniprot_sprot')
        url_unzip(url_sprot, sprot_target)
        log_table[basename(sprot_target)] = date
    except ContentTooShortError:
        print('failed to install {0!s}'.format(basename(sprot_target)))
    try:
        if(not os.path.isdir(PATH_DATABASES+'/uniref90')):
            os.mkdirs(PATH_DATABASES+'/uniref90')
        url_unzip(url_uniref90, uniref90_target)
        log_table[basename(uniref90_target)] = date
    except ContentTooShortError:
        print('failed to install {0!s}'.format(basename(uniref90_target)))
    try:
        safe_retrieve(url_go_pathway, go_pathway_target)
        log_table[basename(go_pathway_target)] = date
    except ContentTooShortError:
        print('failed to install {0!s}'.format(basename(go_pathway_target)))
    try:
        safe_retrieve(url_swiss_enzyme, swiss_enzyme_target)
        log_table[basename(swiss_enzyme_target)] = date
    except ContentTooShortError:
        print('failed to install {0!s}'.format(basename(swiss_enzyme_target)))
    try:
        safe_retrieve(url_pfam_enzyme, pfam_enzyme_target)
        log_table[basename(pfam_enzyme_target)] = date
    except ContentTooShortError:
        print('failed to install {0!s}'.format(basename(pfam_enzyme_target)))
    try:
        url_unzip(url_id_mapping, id_mapping_target)
        log_table[basename(id_mapping_target)] = date
    except ContentTooShortError:
        print('failed to install {0!s}'.format(basename(sprot_map_target)))
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


def main():
    log_table = read_log()
    log_table = download_databases(log_table)
    log_table = subset_dat(id_mapping_target, idmapping_keys, log_table)
    write_log(log_table)


if(__name__ == '__main__'):
    write_log({})
    main()


