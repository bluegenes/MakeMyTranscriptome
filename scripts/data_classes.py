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

    def __init__(self, url, download_path, file_type='', setup_task=None, data_path=None):
        assert file_type in Datafile.valid_file_types
        self.url = url
        self.download_path = download_path
        self.file_type = file_type
        self.setup_task = setup_task
        self.data_file_path = data_path if(data_path is not None) else download_path

    def install(self, force=False):
        if(os.path.exists(self.download_path) and not force):
            return
        else:
            if(self.file_type == ''):
                urlretrieve(self.url, self.download_path)
            elif(self.file_type == 'gz'):
                url_unzip(self.url, self.download_path)
            elif(self.file_type == 'tar'):
                tar_retrieve(self.url, self.download_path)
            else:
                raise Exception('This should never be raised')

    def is_installed(self):
        if(self.setup_task is not None):
            return all(os.path.exists(t) for t in self.setup_task.targets)
        else:
            return os.path.exists(self.download_path)


if(__name__ == '__main__'):
    url_sprot = ('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release'
                 '/knowledgebase/complete/uniprot_sprot.fasta.gz')
    sprot_target = 'uniprot_sprot.fasta'
    sprot_data = Datafile(url_sprot, sprot_target, file_type='gz')
    sprot_data.install(force=False)
