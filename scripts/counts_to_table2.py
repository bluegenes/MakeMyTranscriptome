#Tessa Pierce
#7.9.2015

import sys, re, os

# Input arguments via optparse (works with Python 2.5, which is on the server. Python 2.7 plus prefer argparse...)
from optparse import OptionParser

desc = """ Single Population: many samples --> one counts table"""
parser = OptionParser(description = desc)

### Output Base ###
parser.add_option("--out", "--OutCountsTable", help = "name of the Output Counts Table File" , action="store", type="string", dest="out")

# Directory for Counts files and Output file #
parser.add_option("--inDir", "--INDirectory", help = "path to all input files" , action="store", type="string", dest="inDir", default='./')
parser.add_option("--outDir", "--OUTDirectory", help = "path to all output files" , action="store", type="string", dest="outDir", default='./')
#Threshold #
parser.add_option("-t", "--threshold", type="int", action="store", help = "threshold count number", dest="threshold", default = 2)
#counts and fasta input
parser.add_option("-c", "--counts", "--inCounts", help = "name of the Input counts files" , action="append", type="string", dest="counts", default = [])
#parser.add_option("-f", "--fasta", "--fastaReference", help = "reference fasta file" , action="store", type="string", dest="fasta")
parser.add_option("-g", "--gene-trans-map", "--geneTransMap", help = "reference gene to transcript map" , action="store", type="string", dest="geneTransMap")

# rsem, eXpress or salmon counts files instead of bedtools
parser.add_option("-e", "--eXpress",action="store_true", help = "add this option if the files are from eXpress (false=default)", dest="eXpress", default=False)
parser.add_option("-s", "--salmon",action="store_true", help = "add this option if the files are from salmon (false=default)", dest="salmon", default=False)
parser.add_option("-k", "--kallisto",action="store_true", help = "add this option if the files are from kallisto (false=default)", dest="kallisto", default=False)
parser.add_option("-r", "--rsem",action="store_true", help = "add this option if the files are from rsem (false=default)", dest="rsem", default=False)

(opts, args) = parser.parse_args()
###############################
def countsFromBed(countFile, index, countDt):
    with open(countFile, 'r') as f:
        hitDt = {}
        prevRead = ''
        for line in f:
	    line = line.strip().split('\t')
            hitName = line[0] # contig of the hit
            readName = line[3][:-1] # take off the /1 or /2 from the read name (paired hits get counted as single hit)
            if hitName in hitDt:
                hitDt.get(hitName).update(set([readName]))  # since using a set, pairs will only be represented once.
            else:
                hitDt[hitName] = set([readName])
        f.close()
    for key, val in hitDt.items():
	prevCounts = countDt.get(key)
	prevCounts[index] = str(len(val))
	countDt[key] = prevCounts
    return countDt	

def countsFromExpress(countFile,index, countDt):
    with open(countFile, 'r') as f:
	next(f) #has header 
        for line in f: 
	    line = line.strip().split('\t')	
	    contig = line[1]
	    entry = countDt.get(contig)
	    entry[index] = line[7] # EFFECTIVE counts are index 7 in express file
            #entry[index] = line[6] # ESTIMATED counts are index 6 in express file            
            countDt[contig] = entry
	f.close()
    return countDt

def countsFromSalmon(countFile, index, countDt):
    with open(countFile, 'r') as f:
        for line in f: 
#  	    import pdb;pdb.set_trace()
	    if not line.startswith('#'):
	        line = line.strip().split('\t')	
	        contig = line[0]
	        entry = countDt.get(contig)
                entry[index] = line[3] # want to use NumReads column for input to DESeq2
#		entry[index] = line[2] # if want TPM output
                countDt[contig] = entry
	f.close()
    return countDt

def countsFromKallisto(countFile, index, countDt):
    with open(countFile, 'r') as f:
	next(f) #has header 
        for line in f: 
	    line = line.strip().split('\t')	
	    contig = line[0]
	    entry = countDt.get(contig)
            entry[index] = line[3] # want to use est reads column for input to DESeq2
#	    entry[index] = line[4] # if want TPM output
            countDt[contig] = entry
	f.close()
    return countDt

def countsFromRSEM(countFile, index, countDt):
    with open(countFile, 'r') as f:
	next(f) #has header 
        for line in f: 
	    line = line.strip().split('\t')	
	    contig = line[0] #if isoforms output ... need to check if summing isoforms --> gene file output.
	    entry = countDt.get(contig)
            entry[index] = line[4] # want to use estimated counts column for input to DESeq2
#	    entry[index] = line[5] # if want TPM output
#	    entry[index] = line[6] # if want FPKM output
            countDt[contig] = entry
	f.close()
    return countDt

def getGeneDictionary(transcriptD):
    geneDt = {}
    for key, val in transcriptD.items():
        geneName = val[0]
	if geneName in geneDt:
	    newCounts = [ float(a) + float(b) for a, b in zip(geneDt.get(geneName), val[1:])]
	    geneDt[geneName] = newCounts
	else:
	    geneDt[geneName] = val[1:]
    return geneDt

## main ##
os.chdir(opts.inDir)
countFiles = opts.counts

contigD = {trans: [gene]+ [0]*len(countFiles) for (gene,trans) in (line.strip().split("\t") for line in (open(opts.geneTransMap, 'r')))}

i = 1 # now that gene is the 0th index of the dict, indexing from counts starts at 1 

for file in countFiles:
    if opts.eXpress:
        cD = countsFromExpress(file, i, contigD)
        countFileNames = [name.split('/')[-2] for name in countFiles] #if haven't pulled results.xprs out of it's well-named directory
        #countFileNames = [name.split('/')[-1] for name in countFiles]
    elif opts.salmon:
    	cD = countsFromSalmon(file, i, contigD)
        countFileNames = [name.split('/')[-2] for name in countFiles] #if haven't pulled quant.sf out of it's well-named directory
    elif opts.kallisto:
    	cD = countsFromKallisto(file, i, contigD)
        countFileNames = [name.split('/')[-2] for name in countFiles] #if haven't pulled abundance.txt out of it's well-named directory
    elif opts.rsem:
    	cD = countsFromRSEM(file, i, contigD)
        countFileNames = [name.split('/')[-2] for name in countFiles] #if haven't pulled samples.isoforms.results  out of it's well-named directory
    else:
        cD = countsFromBed(file, i, contigD)
	countFileNames = [name.split('/')[-1] for name in countFiles]
    i = i + 1

geneD = getGeneDictionary(cD)

####### output #######
outD = opts.outDir
if not os.path.exists(outD):
    os.makedirs(outD)

os.chdir(outD)
if opts.eXpress:
    #outName = opts.out + "_estCounts"
    outName = opts.out
else:
    outName = opts.out

outCountTable = open(outName + '.countsTable', 'w')
outThreshTable = open(outName + '_threshold.countsTable', 'w')
outGeneTable = open(outName + '_byGene.countsTable', 'w')

#print header list
outCountTable.write('Contig' + '\t' 'Gene' + '\t' + '\t'.join(countFileNames) + '\n')
outThreshTable.write('Contig'+ '\t' 'Gene' + '\t' + '\t'.join(countFileNames) + '\n')
outGeneTable.write('Contig'  + '\t' + '\t'.join(countFileNames) + '\n')

#write files
# also output a GENE version? Collapse transcripts --> genes!  

for key in sorted(cD): 
    val = contigD.get(key) 
    outCountTable.write(key + '\t' + '\t'.join(map(str,val)) + '\n')
    if all(item >= int(opts.threshold) for item in map(float,val[1:])): # [1:] to take off the gene name
        outThreshTable.write(key + '\t' + '\t'.join(map(str,val[1:])) + '\n') #take off gene name

for key in sorted(geneD):
    val = geneD.get(key)
    outGeneTable.write(key + '\t' + '\t'.join(map(str,val)) + '\n')

outCountTable.close()
outThreshTable.close()
outGeneTable.close()



