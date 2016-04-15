import os
import sys
import gzip
import tarfile
import argparse
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
    tfile.close()
    os.remove(target+'.tar.gz')


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser(description=(
        'Script that allows for the downloading files with the option of unzipping them.'))
    parser.add_argument('url', help=(
        'The address of the file that needs to be downloaded.'))
    parser.add_argument('-t', '--target', default=None, help=(
        'the location to install the file to. Defaults to basename(url) with '
        'compression type trimmed.'))
    parser.add_argument('--type', default='.', help=(
        'the type of decompression that should be performed, if any. '
        'Currently supported types are ".gz" and ".tar.gz".'))
    args = parser.parse_args()
    if(args.target is None):
        args.target = os.path.basename(args.url)
        if(args.target.endswith(args.type)):
            args.target = args.target[:-1*len(args.type)]
    if(args.type == '.'):
        urlretrieve(args.url, args.target)
    elif(args.type == '.gz'):
        url_unzip(args.url, args.target)
    elif(args.type == '.tar.gz'):
        tar_retrieve(args.url, args.target)
    else:
        raise Exception('Unrecognized --type argument : {}'.format(args.type))
