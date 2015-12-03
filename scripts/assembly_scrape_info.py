##############################################
# Tessa Pierce (bluegenes)
# 11.19.2015
# file to get summary statistics out of the assembly quality assessment
##############################################

#edits to make --> use os for all dir path manipulations
#
import os, glob, re, argparse
from parse_qual_metrics import *
import json # best way to write json?
#import assembly_stats --> can just use transrate output for this now.

parser = argparse.ArgumentParser(description="Description: Gather stats on assemblies performed by MMT.")
#required
parser.add_argument('-a','--assemblies_home', help='The directory that contains the assembly directories.')
parser.add_argument('-o','--outBase', help='output assembly stats basename. Results will print to outBase.csv and outBase.json')
#optional paths --> to allow potential directory structure changes
#change this default after the re-factor to make "assess/qualty" it's own tool
parser.add_argument('-q','--qualityDirName', help='Name of directory, within each assembly directory, that contains the quality assessment information.', default='assembly_files')
parser.add_argument('-t','--transrateDirName', help='name of transrate dir', default='transrate_output')
parser.add_argument('-b','--buscoDirName', help='name of transrate dir', default='run_busco_metazoa')
args = parser.parse_args()
assembDir = args.assemblies_home 
outN = args.outBase

# instead of arguments, above, just set the directories as variables --> they should always be in whatever structure *we* determine.
fastqc_pre_trim = os.path.join(assemblies_home, outN, assembly_files,fastqc_pre_trimming)
#fastqc_post_trim = os.path.join(assemblies_home, outN, assembly_files,fastqc_pre_trimming)
fastqc_final = os.path.join(assemblies_home, outN, assembly_files,fastqc_final_reads_paired)

# busco can be run with any of 5 (6?) databases --> need to handle all naming options... split name on last '_', grab database name for the dictionary.


def is_assembly_dir(assembliesHome, d, quality_dir, transrate_dir):
    '''Check that at least a fasta file and a transrate report are available in the directory'''
    assembPath = assembliesHome + '/' + d 
    if (len(glob.glob(assembPath + '/*.fasta')) > 0 and len(glob.glob(assembPath + '/' + quality_dir + '/' + transrate_dir + '/assemblies.csv')))> 0:
        return True
    else:
        return False

def scrape(assembliesHome, d ,quality_dir, transrate_dir, busco_dir):
    '''	Extracts information from all of the relevant files in an assembly directory.
    IF ADDING ANOTHER ASSESSMENT TOOL: Add a line for the tool here; add a file parser to parse_qual_metrics.py '''
    assemblyPath = assembliesHome + '/' + d + '/' + quality_dir
    infoD = {}
    infoD.update(transrate_parser(glob.glob(assemblyPath + '/' + transrate_dir +"/*assemblies.csv")))
    infoD.update(busco_parser(glob.glob(assemblyPath + "/" + busco_dir + "/short_summary*")))
    infoD.update(transrate_readcount_parser(glob.glob(assemblyPath + "/" + transrate_dir + '/' + d + '/' + "*read_count.txt")))
    infoD.update(cegma_parser(glob.glob(assemblyPath +"/*.completeness_report")))
    infoD.update(detonate_parser(glob.glob(assemblyPath +"/*.score")))
    infoD.update(fastqc_parser(glob.glob(assemblyPath + fastqc_dir + "/fastqc_data.txt")))
    #don't think we need these anymore ... transrate gets most of the stats for us
    #infoD.update(fasta_parser(glob.glob(assemblyPath + "/*.fasta")))
    #infoD.update(stats_parser(glob.glob(assemblyPath +"/*.stats")))
    return infoD

def writeAssemblyStats_csv(statsDict, outFi):
    outF = open(outFi + '.csv', 'w')
    outS = open(outFi + '_subset.csv', 'w')
    header = True
    for assemblyName in sorted(statsDict):
    	stats = statsDict[assemblyName]
    	statsKeys = list(stats.keys())
	statsVals = list(stats.values())
	subsetKeys= ['readCount', 'n_seqs', 'n50','metazoa_%_complete', 'n_bases', 'gc', 'p_good_mapping']
	if header:
	    outF.write("assemblyName" + ',' + ','.join(map(str,statsKeys)) + '\n') 
	    outS.write("assemblyName" + ',' + ','.join(map(str,subsetKeys)) + '\n')
	    header=False
	outF.write(assemblyName + ',' + ','.join(map(str,statsVals)) + '\n')#giant tsv for now: 
	# make a much smaller abbreviated file with *just* the params I want to graph (for now). should be able to pull out specific params from dict
	subsetD = [stats[x] for x in subsetKeys]
	outS.write(assemblyName + ',' + ','.join(map(str,subsetD)) + '\n')

def writeAssemblyStats_json(statsDict, outF):
    outF = open(outF + '.json', 'w')
    # write out ...

# main #
outFile = assembDir + '/' + outN  #general output filename (for both .csv and .json) 
assembliesD = {} 

for d in os.walk(assembDir).next()[1]: # list *only* directories
    if is_assembly_dir(assembDir,d, args.qualityDirName, args.transrateDirName):
        assembliesD[d] = scrape(assembDir, d, args.qualityDirName, args.transrateDirName, args.buscoDirName)
	print assembliesD
	#import pdb;pdb.set_trace()
    writeAssemblyStats_csv(assembliesD, outFile)
#    writeAssemblyStats_json(assembliesD, outFile)


#nolan's write info:
#        for key in sorted(results[assembly]):
# 	    print('\t'+str(key)+' : '+str(results[assembly][key]))




