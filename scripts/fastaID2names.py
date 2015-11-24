#!/matta1/biotools/anaconda/bin/python

# script to get id: full description from the fasta file (swissprot, nr, ...)
# update: now gets length of fasta entry as well. Note, this is *protein* length (# amino acids)

#simple version, doesn't use fasta parser ---> potentially edit to use Nolan's fasta parser

import sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', help='tab-separated database lookup: full name file for reference (eg nr or swissprot)')
args = parser.parse_args()


outF = open(args.fasta.split('.fa')[0] + '.id2names', 'w')

bp = 0
name =''
with open(args.fasta) as f:
    for line in f:
        if line.startswith(">"):
	    if len(name) > 0: # need to write *after* get bp len
                outF.write( ID + "\t" + name + '\t' + str(bp) + '\n')
            	bp = 0
	    entry = line.rstrip().split(" ", 1)
            ID = entry[0][1:]
	    name = entry[1]
        else:
            bp = bp + len(line.rstrip())
    #catch last contig
    outF.write( ID + "\t" + name + '\t' + str(bp) + '\n')

    outF.close()
    f.close()



