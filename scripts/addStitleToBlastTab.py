import sys, argparse

parser = argparse.ArgumentParser()

parser.add_argument('--db2Name', help='tab-separated database lookup: full name file for reference (eg nr or swissprot)')
parser.add_argument('-b','--blast', help='blast input file')

args = parser.parse_args()

blastFile = []
hitList = []
hitDt = {}
#header = ['query_id', 'subject_id', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bitscore', 'full_name', 'subject_length']
header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle', 'slen']


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

sys.stdout.write('\t'.join(header) + '\n')

for line in blastFile:
    extraInfo = hitDt.get(line[1], [])
    line = line + extraInfo
    sys.stdout.write('\t'.join(map(str,line)) + '\n')

