#!/matta1/biotools/anaconda/bin/python


# script to get id: full description from the fasta file (swissprot, nr, ...)
# usage: grep ">" uniprot_sprot.fasta | fastaID2names.py > sprot_id2Name.txt


import sys

for line in sys.stdin:
    if line.startswith(">"):
        entry = line.rstrip().split(" ", 1)
        ID = entry[0][1:]
	    name = entry[1]
        sys.stdout.write( ID + "\t" + name + '\n')

sys.stdin.close()
sys.stdout.close()
sys.exit()


