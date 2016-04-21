#!/matta1/biotools/anaconda/bin/python

# script to get id: full description from the fasta file (swissprot, nr, ...) + length of fasta entry as well. Note, this is *protein* length (# amino acids)

## this one needs to read in a gz file, so we don't have to unzip it.

import sys, argparse, gzip

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', help='tab-separated database lookup: full name file for reference (eg nr or swissprot)')
args = parser.parse_args()


#outF = open(args.fasta.split('.fa')[0] + '.id2names', 'w')
outF = open(args.fasta.split('.gz')[0] + '.stitle', 'w')

bp = 0
name =''
with gzip.open(args.fasta) as f:
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



