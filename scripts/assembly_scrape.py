##############################################
# Tessa Pierce (bluegenes)
# file to get summary statistics out of the assembly quality assessment
##############################################

#
import os, glob, re, argparse
from parse_qual_metrics import *
import json # best way to write json?
import pandas as pd
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import seaborn as sns

#import assembly_stats --> can just use transrate output for this now.

parser = argparse.ArgumentParser(description="Description: Gather stats on assemblies performed by MMT.")
parser.add_argument('-a','--assemblies_home', help='The directory that contains the assembly directories.')
parser.add_argument('-o','--outBase', help='output assembly stats basename. Results will print to outBase.csv and outBase.json')

parser.add_argument('-q','--qualityDir', help='optional: alternative name for quality directory', default='quality_files')
parser.add_argument('-t','--transrateDir', help='optional: alternative name for transrate directory', default='transrate_quality')
parser.add_argument('--transratePostAssemblyDir', help='optional: alternative name for transrate directory', default='transrate_post_assembly')

args = parser.parse_args()
assembDir = args.assemblies_home 
outN = args.outBase
###########################################
# set up our directory stucture, filenames
###########################################
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

fileD = {'cegma':cegmaFile, 'busco':buscoFile, 'transrate':transrateFile, 'transrateRC': transrateReadCountF, 'detonate':detonateFile,'fastqc':fastqcFile}
###########################################

def is_assembly_dir(assembliesHome, d, quality_dir, transrate_dir):
    '''Check that at least a fasta file and a transrate report are available in the directory'''
    fastaF = os.path.join(assembliesHome, d, '*.fasta')
    print fastaF
    transrateF = os.path.join(assembliesHome, d, quality_dir, transrate_dir, transrateFile)
    print transrateF
    if (len(glob.glob(fastaF)) > 0 and len(glob.glob(transrateF)))> 0:
        print fastaF
        print d
        return True
    else:
        return False

def writeAssemblyStats_csv(statsDF, outFi):
    outF = open(outFi + '.csv', 'w') #full
    statsDF.to_csv(outF, sep=',')
   #outS = open(outFi + '_subset.csv', 'w') #subset of info
#    header = True
#    for assemblyName in sorted(statsDict):
#    	stats = statsDict[assemblyName]
#    	statsKeys = list(stats.keys())
#	statsVals = list(stats.values())
	#subsetKeys= ['Transrate_readCount', 'n_seqs', 'n50','metazoa_%_complete', 'n_bases', 'gc', 'p_good_mapping']
#	if header:
#	    outF.write("assemblyName" + ',' + ','.join(map(str,statsKeys)) + '\n') 
	   # outS.write("assemblyName" + ',' + ','.join(map(str,subsetKeys)) + '\n')
#	    header=False
#	outF.write(assemblyName + ',' + ','.join(map(str,statsVals)) + '\n')#giant tsv for now: 
	# make a much smaller abbreviated file with *just* the params I want to graph (for now). should be able to pull out specific params from dict
	#subsetD = [stats[x] for x in subsetKeys]
	#outS.write(assemblyName + ',' + ','.join(map(str,subsetD)) + '\n')
#def writeAssemblyStats_json(statsDict, outF):
#    outF = open(outF + '.json', 'w')
    # write out ...

# main #
outFile = os.path.join(assembDir, outN)  #general output filename (for both .csv and .json) 
assembliesD = {} 

def scrape_single_assembly(args, assembPath, d, fileD):
    ''' Extracts information from all of the relevant files in an assembly directory.
    IF ADDING ANOTHER ASSESSMENT TOOL: Add a line for the tool here; add a file parser to parse_qual_metrics.py '''
    #set up dir paths
    qualityPath = os.path.join(assembPath, args.qualityDir)
    transratePath = os.path.join(qualityPath, args.transrateDir)
    #scrape info into a dict
    infoD = {}
    for prog, fileName in fileD.items():
        if prog == 'transrate':
    	    transrateD = transrate_parser(transratePath, fileName) # the d from here can be diff btwn assemblies, depending on whether reads were used or not...
	    infoD.update(transrateD)
	if prog == 'busco':
    	    buscoDirs = glob.glob(os.path.join(qualityPath, 'run_*'))
    	    buscoD = {}
    	    for buscoDir in buscoDirs:
       		buscoD.update(busco_parser(os.path.join(qualityPath, buscoDir), fileName))
	    infoD.update(buscoD)
	if prog == 'cegma':
	    cegmaD = cegma_parser(qualityPath, fileName)
	    infoD.update(cegmaD)
	if prog == 'detonate':
	    detonateD = detonate_parser(qualityPath, fileName)
	    infoD.update(detonateD)
	if prog == 'fastqc':
	    fastqcPreTrimD = fastqc_parser(os.path.join(assemblyPath, fastqc_pre_trim), fileName)
	    fastqcFinalD = fastqc_parser(os.path.join(assemblyPath, fastqc_final), fileName)
	    infoD.update(fastqcPreTrimD)
	    infoD.update(fastqcFinalD)
    print d
    print infoD
    return infoD

    
#in progress
#def plot_assembly_comparisons(assembDF):
#some plotting defaults
#    sns.set(style="white")
#    sns.despine()
#    assembOrder = ['garlic', 'garlic_trunc200','garlic_trunc150','garlic_trunc100','garlic_trunc50']
#    meltedDF = pd.melt(assembDF)
#    g = sns.factorplot(x="assembly", y="value", col="variable",ci=None, x_order=assembOrder, col_wrap=2, kind="bar", data=meltedDF, size=3, aspect=2, sharex=False,sharey=False)
#    g = sns.factorplot(x="assembly", y="value", col="variable",ci=None, col_wrap=2, kind="bar", data=meltedDF, size=3, aspect=2, sharex=False,sharey=False)
#    g.savefig('transrate_plots.png')




assembliesD = {}
for d in os.walk(assembDir).next()[1]: # list *only* directories
    if is_assembly_dir(assembDir,d, args.qualityDir, args.transrateDir):
        assemblyPath = os.path.join(assembDir, d)
	assemblyName = d # just being explicit about where we're grabbing the name from each assembly
	assembliesD[assemblyName] = scrape_single_assembly(args, assemblyPath, d, fileD) 
    print assembliesD
#assembliesDF = combineDtoDF(assembliesD)
assembliesDF = pd.DataFrame.from_dict(assembliesD)
#    plot_assembly_comparisons(assembliesDF)
writeAssemblyStats_csv(assembliesDF, outFile)
    #writeAssemblyStats_json(assembliesDF, outFile)






	
''' unused 
#	transrateD =transrate_parser(transratePath, transrateFile)
#	infoD.update(transrate_parser(transratePath, transrateFile))
#	buscoDirs = glob.glob(os.path.join(qualityPath, 'run_*'))
	#for buscoD in buscoDirs:
   	#    infoD.update(busco_parser(os.path.join(qualityPath, buscoD), buscoFile))
#	infoD.update(transrate_readcount_parser(os.path.join(qualityPath, transrateDirName, d), transrateReadCountF))
#        infoD.update(cegma_parser(qualityPath, cegmaFile))
#        infoD.update(detonate_parser(qualityPath, detonateFile))
#        infoD.update(fastqc_parser(os.path.join(assemblyPath, fastqc_pre_trim), fastqcFile))
#        infoD.update(fastqc_parser(os.path.join(assemblyPath, fastqc_final), fastqcFile))

#don't actually need this unless dicts within dicts within dicts...
#def combineDtoDF(assembliesD):
   #make sense of it all. Avoid any problems with transrateD's being different, etc.
    #assembliesDF ={}
#    names = []
#    infos = []
#    for aName, dt in assembliesD.iteritems():
    
#        names.append(aName)
#        infos.append(pd.DataFrame.from_dict(dt, orient='index'))
#        print names
#        print infos
#    assembliesDF = pd.concat(infos, keys=names)
#    print assembliesDF
#    return assembliesDF

'''
