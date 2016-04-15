import os
import mmt_defaults
import json_config

'''
import os
import sys
import gzip
import tarfile
if(sys.version[0] == '3'):
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve


def url_unzip(source, target):
    urlretrieve(source, target+'.gz')
    f = gzip.open(target+'.gz', 'rb')
    g = open(target, 'wb')
    for line in f:
        g.write(line)
    f.close()
    g.close()
    os.remove(target+'.gz')


def tar_retrieve(source, target):
    urlretrieve(source, target+'.tar.gz')
    tfile = tarfile.open(target+'.tar.gz', 'r:gz')
    tfile.extractall(target)
    os.remove(target+'.tar.gz')


class Datafile:

    valid_file_types = set(['', 'gz', 'tar'])

    def __init__(self, url, download_path, file_type=''):
        assert file_type in Datafile.valid_file_types
        self.url = url
        self.download_path = download_path
        self.file_type = file_type
        self.db_path = download_path

    def install(self, force=False):
        if(force or not self.is_installed()):
            if(self.file_type == ''):
                urlretrieve(self.url, self.download_path)
            elif(self.file_type == 'gz'):
                url_unzip(self.url, self.download_path)
            elif(self.file_type == 'tar'):
                tar_retrieve(self.url, self.download_path)
            else:
                raise Exception('This should never be raised')

    def is_installed(self):
        return os.path.exists(self.download_path)


class Database:

    def __init__(self, df, setup_task, db_path):
        self.datafile = df
        self.setup_task = setup_task
        self.db_path = db_path

    def install(self, force=False):
        self.datafile.install(force)
        if(force or not self.is_installed()):
            self.setup_task.run()

    def is_installed(self):
        return all(os.path.exists(t) for t in self.setup_task.targets)


url_sprot = ('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release'
             '/knowledgebase/complete/uniprot_sprot.fasta.gz')
sprot_target = 'uniprot_sprot.fasta'
sprot_data = Datafile(url_sprot, sprot_target, file_type='gz')
sprot_data.install(force=False)
'''


def file_exists(obj):
    if(os.path.exists(obj)):
        return obj
    else:
        raise Exception('Unable to locate {0!s}.'.format(obj))


def get_db_config():
    db_config = json_config.json_config()
    db_config.add_config('uniprot_sprot_fasta',
                         type=file_exists,
                         default=mmt_defaults.PATH_UNIPROT_SPROT)
    db_config.add_config('uniref90_fasta',
                         type=file_exists,
                         default=mmt_defaults.PATH_UNIREF90)
    db_config.add_config('nr_fasta',
                         type=file_exists,
                         default=mmt_defaults.PATH_NR)
    db_config.add_config('busco_metazoa',
                         type=file_exists,
                         default=mmt_defaults.PATH_BUSCO_METAZOA)
    db_config.add_config('busco_anthropoda',
                         type=file_exists,
                         default=mmt_defaults.PATH_BUSCO_ANTHROPODA)
    db_config.add_config('busco_vertebrata',
                         type=file_exists,
                         default=mmt_defaults.PATH_BUSCO_VERTEBRATA)
    db_config.add_config('busco_eukaryota',
                         type=file_exists,
                         default=mmt_defaults.PATH_BUSCO_EUKARYOTA)
    db_config.add_config('busco_fungi',
                         type=file_exists,
                         default=mmt_defaults.PATH_BUSCO_FUNGI)
    db_config.add_config('busco_bacteria',
                         type=file_exists,
                         default=mmt_defaults.PATH_BUSCO_BACTERIA)
    db_config.add_config('busco_plant',
                         type=file_exists,
                         default=mmt_defaults.PATH_BUSCO_PLANT)
    db_config.add_config('pfam',
                         type=file_exists,
                         default=mmt_defaults.PATH_PFAM_DATABASE)
    db_config.add_config('go_pathway',
                         type=file_exists,
                         default=mmt_defaults.PATH_GO_PATHWAY)
    db_config.add_config('id_mapping',
                         type=file_exists,
                         default=mmt_defaults.PATH_ID_MAPPING)
    db_config.add_config('uniprot_sprot_map',
                         type=file_exists,
                         default=mmt_defaults.PATH_UNIPROT_SPROT_MAP)
    db_config.add_config('kog_functional',
                         type=file_exists,
                         default=mmt_defaults.PATH_KOG_FUNCTIONAL)
    db_config.add_config('goslim_generic',
                         type=file_exists,
                         default=mmt_defaults.PATH_SLIM_GENERIC)
    db_config.add_config('idmapping_selected',
                         type=file_exists,
                         default=mmt_defaults.PATH_IDMAPPING_SELECTED)
    db_config.add_config('pfam_enzyme',
                         type=file_exists,
                         default=mmt_defaults.PATH_PFAM_ENZYME)
    return db_config


class database:

    def __init__(self, url, download_location, call_path=None, ftype='.'):
        self.url = url
        self.download_location = download_location
        self.call_path = call_path if(call_path is not None) else download_location
        self.type = ftype


def get_dbs(json_config=mmt_defaults.PATH_DB_CONFIG_FILE, defaults=False):
    ret = {}
    db_conf = get_db_config()
    config = db_conf.defaults if(defaults) else db_conf.load_config(json_config, True)

    ret['uniprot_sprot'] = database(
        mmt_defaults.URL_UNIPROT_SPROT,
        config['uniprot_sprot_fasta'],
        os.path.join(mmt_defaults.PATH_UNIPROT_SPROT_DIR,
                     os.path.basename(config['uniprot_sprot_fasta'])),
        '.gz')

    ret['uniref90'] = database(
        mmt_defaults.URL_UNIREF90,
        config['uniref90_fasta'],
        os.path.join(mmt_defaults.PATH_UNIREF90_DIR,
                     os.path.basename(config['uniref90_fasta'])),
        '.gz')

    ret['nr'] = database(
        mmt_defaults.URL_NR,
        config['nr_fasta'],
        os.path.join(mmt_defaults.PATH_NR_DIR,
                     os.path.basename(config['nr_fasta'])),
        '.gz')

    ret['busco_metazoa'] = database(
        mmt_defaults.URL_BUSCO_METAZOA,
        config['busco_metazoa'],
        '.tar.gz')

    ret['busco_anthropoda'] = database(
        mmt_defaults.URL_BUSCO_ANTHROPODA,
        config['busco_anthropoda'],
        '.tar.gz')

    ret['busco_vertebrata'] = database(
        mmt_defaults.URL_BUSCO_VERTEBRATA,
        config['busco_vertebrata'],
        '.tar.gz')

    ret['busco_eukaryota'] = database(
        mmt_defaults.URL_BUSCO_EUKARYOTA,
        config['busco_eukaryota'],
        '.tar.gz')

    ret['busco_fungi'] = database(
        mmt_defaults.URL_BUSCO_FUNGI,
        config['busco_fungi'],
        '.tar.gz')

    ret['busco_bacteria'] = database(
        mmt_defaults.URL_BUSCO_BACTERIA,
        config['busco_bacteria'],
        '.tar.gz')

    ret['busco_plant'] = database(
        mmt_defaults.URL_BUSCO_PLANT,
        config['busco_plant'],
        '.tar.gz')

    ret['go_pathway'] = database(
        mmt_defaults.URL_GO_PATHWAY,
        config['go_pathway'])

    ret['pfam_enzyme'] = database(
        mmt_defaults.URL_PFAM_ENZYME,
        config['pfam_enzyme'])

    ret['id_mapping'] = database(
        mmt_defaults.URL_ID_MAPPING,
        config['id_mapping'])

    ret['uniprot_sprot_map'] = database(
        mmt_defaults.URL_UNIPROT_SPROT_MAP,
        config['uniprot_sprot_map'])

    ret['kog_functional'] = database(
        mmt_defaults.URL_KOG_FUNCTIONAL,
        config['kog_functional'])

    ret['goslim_generic'] = database(
        mmt_defaults.URL_SLIM_GENERIC,
        config['goslim_generic'])

    ret['pfam'] = database(
        mmt_defaults.URL_PFAM_DATABASE,
        config['pfam'],
        os.path.join(mmt_defaults.PATH_PFAM_DIR, os.path.basename(config['pfam'])))

    ret['idmapping_selected'] = database(
        mmt_defaults.URL_IDMAPPING_SELECTED,
        config['idmapping_selected'])

    return ret


if(__name__ == '__main__'):
    pass
