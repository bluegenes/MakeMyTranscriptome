# author: bluegenes
import os
import argparse
import sys
import warnings


def quants_to_dict(quantFiles, tpm_index):
    tpmD = {}
    for q in quantFiles:
        with open(q, 'r') as f:
            for line in f:
                if (line.startswith('Name') or line.startswith('#')):
                    continue
                line = line.strip().split('\t')
                contig = line[0]
                prev_tpm = tpmD.get(contig, 0)
                tpm = prev_tpm + float(line[tpm_index])
                tpmD[contig] = tpm
    return tpmD


def main(assembly, quant_files, tpm_threshold, out, tpm_col_index):
    if len(quant_files) == 0:
        warnings.warn('No quant files passed in; '+ assembly + 'cannot be filtered.')
    else:
        quantD = quants_to_dict(quant_files, tpm_col_index)
        filteredSet = set([contig for contig,tpm in quantD.items() if tpm >= tpm_threshold])
        with open(assembly, 'r') as f:
            with open(out, 'w') as o:
	        for line in f:
                    if line.startswith('>'):
                        write_bases = False
                        contig = line.rstrip()[1:].split(' ')[0] # salmon names = only before the first space 
                        if contig in filteredSet:
                            o.write(line)
                            write_bases=True
                    else:
                        if write_bases:
                            o.write(line)


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument('--assembly', action='store', type=str)
    parser.add_argument('--tpm', action='store', type=float,default=0.5)
    parser.add_argument('--quant_files', action='append', default=[])
    parser.add_argument('-o', '--out', action='store', type=str)
    parser.add_argument('--tpm_column_index', action='store', type=int,default=3)
    args = parser.parse_args()
    main(args.assembly, args.quant_files, args.tpm, args.out, args.tpm_column_index)

