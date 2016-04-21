import argparse
import os


def subset_dat(dat_file, biocyc, eggnog, ko, orthodb):
    idmapping_keys = {'BioCyc': biocyc,
                      'eggNOG': eggnog,
                      'KO': ko,
                      'OrthoDB': orthodb}
    dat = open(dat_file)
    key_file = {key: open(idmapping_keys[key], 'w') for key in idmapping_keys}
    for line in dat:
        values = line.split('\t')
        if(values[1] in key_file):
            key_file[values[1]].write(line)
    for key in key_file:
        key_file[key].close()


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument('dat_file', help='the dat file we wish to subset')
    parser.add_argument('--biocyc',
                        help='the path where biocyc entries should be written to.',
                        default=os.devnull)
    parser.add_argument('--eggNOG',
                        help='the path where eggnog entries should be written to.',
                        default=os.devnull)
    parser.add_argument('--ko',
                        help='the path where eggnog entries should be written to.',
                        default=os.devnull)
    parser.add_argument('--orthodb',
                        help='the path where orthodb entries should be written to',
                        default=os.devnull)
    args = parser.parse_args()
    subset_dat(args.dat_file, args.biocyc, args.eggnog, args.ko, args.orthodb)
