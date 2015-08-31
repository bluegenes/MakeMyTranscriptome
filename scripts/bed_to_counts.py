#Tessa Pierce
#6.2.2103

#Parse BED file into a Tab-separated file containing contig name and count data

# 0 = contig name
# 1 = hit start pos
# 2 = hit end pos
# 3 = read name
# 4
# 5
# 6 = contig name
# 7 = contig start pos
# 8 = contig end pos


# at this point, we don't care about the fraction of overlap -- any "hit" counts. so, the important features are the contig name (0) and the read name (3). Since we have some paired end alignments and some discordantly mapped pairs (~single end alignments), we can't use the bedtools map pairs (pairtobed) function. However, since the reads of a pair are not independent, we should treat any hit to either read of the pair as one hit, and if they hit the same gene of interest (GOI), count it as only 1 hit.

import sys, re

basename = sys.argv[1]

inBED = open(basename + '.bed', 'r')
outCOUNTS = open(basename + '.counts', 'w')

Bedlines = [ x.strip().split('\t') for x in inBED.readlines() ]


hitDt = {}
prevRead = ''

#line = inBED.readline().strip().split('\t')

for line in Bedlines:     
    hitName = line[0]
    readName = line[3][:-1]
    #contigLen = float(line[8]) - float(line[7])
    contigLen = float(line[7]) - float(line[6])
    if readName != prevRead:
        prevRead = readName
	if hitName in hitDt:
	    counts = hitDt.get(hitName) + 1
        else: 
	    counts = 1
	hitDt[hitName] = counts


for key, val in hitDt.items():
    outCOUNTS.write(str(key) + '\t' + str(val) + '\n' )


inBED.close()
outCOUNTS.close()

