##############################################
# Tessa Pierce (bluegenes)
# file to get summary statistics out of the assembly quality assessment
##############################################

#edits to make --> use os for all dir path manipulations
#
import os, glob, re, argparse
from parse_qual_metrics import *
import json # best way to write json?
#import assembly_stats --> can just use transrate output for this now.

parser = argparse.ArgumentParser(description="Description: Gather stats on assemblies performed by MMT.")
parser.add_argument('-a','--assemblies_home', help='The directory that contains the assembly directories.')
parser.add_argument('-o','--outBase', help='output assembly stats basename. Results will print to outBase.csv and outBase.json')

args = parser.parse_args()
assembDir = args.assemblies_home 
outN = args.outBase
###########################################
# set up our directory stucture, filenames
###########################################
qualityDirName = "quality_files"
transrateDirName = "transrate_output"
buscoDirBase = "run_busco_*"
fastqc_pre_trim = "fastqc_pre*"
fastqc_final = "fastqc_final*"
#third fastqc if subsampling?

cegmaFile = "*.completeness_report"
buscoFile = "short_summary*"
transrateFile = "*assemblies.csv"
transrateReadCountF = "*read_count.txt" # need readcount? 
detonateFile = "*.score"
fastqcFile = "fastqc_data.txt"
###########################################

def is_assembly_dir(assembliesHome, d, quality_dir, transrate_dir):
    '''Check that at least a fasta file and a transrate report are available in the directory'''
    fastaF = os.path.join(assembliesHome, d, '*.fasta')
    transrateF = os.path.join(assembliesHome, d, quality_dir, transrate_dir,'assemblies.csv')
    print transrateF
    if (len(glob.glob(fastaF)) > 0 and len(glob.glob(transrateF)))> 0:
        print fastaF
        print d
        return True
    else:
        return False

def writeAssemblyStats_csv(statsDict, outFi):
    outF = open(outFi + '.csv', 'w') #full
    #outS = open(outFi + '_subset.csv', 'w') #subset of info
    header = True
    for assemblyName in sorted(statsDict):
    	stats = statsDict[assemblyName]
    	statsKeys = list(stats.keys())
	statsVals = list(stats.values())
	#subsetKeys= ['Transrate_readCount', 'n_seqs', 'n50','metazoa_%_complete', 'n_bases', 'gc', 'p_good_mapping']
	if header:
	    outF.write("assemblyName" + ',' + ','.join(map(str,statsKeys)) + '\n') 
	   # outS.write("assemblyName" + ',' + ','.join(map(str,subsetKeys)) + '\n')
	    header=False
	outF.write(assemblyName + ',' + ','.join(map(str,statsVals)) + '\n')#giant tsv for now: 
	# make a much smaller abbreviated file with *just* the params I want to graph (for now). should be able to pull out specific params from dict
	#subsetD = [stats[x] for x in subsetKeys]
	#outS.write(assemblyName + ',' + ','.join(map(str,subsetD)) + '\n')
#def writeAssemblyStats_json(statsDict, outF):
#    outF = open(outF + '.json', 'w')
    # write out ...

# main #
outFile = assembDir + '/' + outN  #general output filename (for both .csv and .json) 
assembliesD = {} 

for d in os.walk(assembDir).next()[1]: # list *only* directories
    if is_assembly_dir(assembDir,d, qualityDirName, transrateDirName):
        assemblyPath = os.path.join(assembDir, d)
        qualityPath = os.path.join(assembDir, d, qualityDirName)
        ''' Extracts information from all of the relevant files in an assembly directory.
        IF ADDING ANOTHER ASSESSMENT TOOL: Add a line for the tool here; add a file parser to parse_qual_metrics.py '''
        infoD = {}
	infoD.update(transrate_parser(os.path.join(qualityPath, transrateDirName), transrateFile))
	buscoDirs = glob.glob(os.path.join(qualityPath, 'run_*'))
	for buscoD in buscoDirs:
   	    infoD.update(busco_parser(os.path.join(qualityPath, buscoD), buscoFile))
        infoD.update(transrate_readcount_parser(os.path.join(qualityPath, transrateDirName, d), transrateReadCountF))
        infoD.update(cegma_parser(qualityPath, cegmaFile))
        infoD.update(detonate_parser(qualityPath, detonateFile))
        infoD.update(fastqc_parser(os.path.join(assemblyPath, fastqc_pre_trim), fastqcFile))
        infoD.update(fastqc_parser(os.path.join(assemblyPath, fastqc_final), fastqcFile))
	assemblyName = d # just being explicit about where we're grabbing the name from each assembly
     	print assemblyName
	print infoD  
	assembliesD[assemblyName] = infoD
    writeAssemblyStats_csv(assembliesD, outFile)
    #writeAssemblyStats_json(assembliesD, outFile)
	
	




