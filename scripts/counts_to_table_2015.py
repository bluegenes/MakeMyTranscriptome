#Tessa Pierce
#1.7.2104

import sys, re, os

# Input arguments via optparse (works with Python 2.5, which is on the server. Python 2.7 plus prefer argparse...)
from optparse import OptionParser

desc = """ Single Population: many samples --> one counts table
"""
parser = OptionParser(description = desc)

### Output Base ###
parser.add_option("--out", "--OutCountsTable", help = "name of the Output Counts Table File" , action="store", type="string", dest="out")

# Directory for Counts files and Output file #
parser.add_option("--inDir", "--INDirectory", help = "path to all input files" , action="store", type="string", dest="inDir")
parser.add_option("--outDir", "--OUTDirectory", help = "path to all output files" , action="store", type="string", dest="outDir")

#Threshold #
parser.add_option("-t", "--threshold", type="int", action="store", help = "threshold count number", dest="threshold", default = 2)

# eXpress counts files instead of bedtools
parser.add_option("-e", "--eXpress",action="store_true", help = "add this option if the files are from eXpress (false=default)", dest="eXpress", default=False)

parser.add_option("-c", "--counts", "--inCounts", help = "name of the Input counts files" , action="append", type="string", dest="counts", default = [])

parser.add_option("-b", "--bed", "--contigList", help = "file, i.e. a bed file OR an eXpress file results.xprs.genes file containing all contig names in the first column" , action="store", type="string", dest="bed")

parser.add_option("-f", "--fasta", "--fastaReference", help = "reference fasta file" , action="store", type="string", dest="fasta")

(opts, args) = parser.parse_args()

#main
os.chdir(opts.inDir)
countFiles = opts.counts


contigD = {}
with open(opts.fasta) as f:
    for line in f:
        if line.startswith('>'):
	    contigD[line[1:]] = [0]*len(countFiles)

if opts.eXpress:
    contigIndex = 1
    countIndex  = 7
    startIndex = 1
    contigD =  {values[1]: [0]*len(countFiles) for values in (x.split('\t') for x in open(countFiles[0], 'r')) } 
    del contigD["target_id"]
else:
    contigIndex = 0
    startIndex = 0
    countIndex = 1
#    bedF = open(opts.bed, 'r')
#    contigD =  {values[0]: [0]*len(countFiles) for values in (x.strip().split('\t') for x in bedF) } 

i = 0
for file in countFiles:
    with open(file, 'r') as f:
        fLines = [x.strip().split('\t') for x in f]
        for line in fLines[startIndex:]:
	    entry = contigD.get(line[contigIndex])
	    entry[i] = line[countIndex]
	    contigD[line[contigIndex]] = entry
    i = i + 1
    f.close()

####### output #######
outD = opts.outDir
if not os.path.exists(outD):
    os.makedirs(outD)

os.chdir(outD)
outCountTable = open(opts.out + '.countsTable', 'w')
outThreshTable = open(opts.out + '_threshold.countsTable', 'w')
#print header list
outCountTable.write('Contig' +'\t' + '\t'.join(countFiles) + '\n')
outThreshTable.write('Contig'+ '\t' + '\t'.join(countFiles) + '\n')

#write files
for key in sorted(contigD): #trying sorted --> need same order
    val = contigD.get(key)
    outCountTable.write(key + '\t' + '\t'.join(map(str,val)) + '\n')
    if all(item >= 2 for item in map(float,val)):
        outThreshTable.write(key + '\t' + '\t'.join(map(str,val)) + '\n') 

outCountTable.close()
outThreshTable.close()
#bedF.close()



