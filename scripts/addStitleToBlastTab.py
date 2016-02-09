import sys, argparse

parser = argparse.ArgumentParser()

parser.add_argument('--db2Name', help='tab-separated database lookup: full name file for reference (eg nr or swissprot)')
parser.add_argument('-b','--blast', help='blast input file')

args = parser.parse_args()

blastFile = []
hitList = []
hitDt = {}
with open(args.blast, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
	blastFile.append(line)
	hitList.append(line[1])

hitSet = set(hitList)
#potentially huge file --> don't want this in memory
with open(args.db2Name, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
        if line[0] in hitSet:
            hitDt[line[0]] = line[1:]
    f.close()

for line in blastFile:
    extraInfo = hitDt.get(line[1], [])
    line = line + extraInfo
    sys.stdout.write('\t'.join(map(str,line)) + '\n')

