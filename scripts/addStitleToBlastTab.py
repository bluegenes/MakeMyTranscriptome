import sys, argparse

parser = argparse.ArgumentParser()

parser.add_argument('--db2Name', help='tab-separated database lookup: full name file for reference (eg nr or swissprot)')
parser.add_argument('-b','--blast', help='blast input file')
parser.add_argument('-o','--out', help='output file')

args = parser.parse_args()

#can't do this bc need to preserve order:
#blastD = {x[2]:x for x in args.blast.rstrip().split('\t')}

blastOrder = []
blastD = {}
with open(args.blast, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
	#import pdb; pdb.set_trace()
	blastOrder.append(line[1])
	blastD[line[1]] = line
    f.close()

#potentially huge file --> don't want this in memory
with open(args.db2Name, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
	hitInfo = blastD.get(line[0], None)
	if hitInfo is not None:
	    hitInfo.append(line[1])
    f.close()

outExtendedTab = open(args.out, 'w')
for hit in blastOrder:
    outExtendedTab.write('\t'.join(map(str,blastD[hit])) + '\n')




