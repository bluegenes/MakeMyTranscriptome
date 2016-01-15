##############################################
# Tessa Pierce
# 1/2016
# scrape assembly quality stats
##############################################

import os, glob, re, argparse
from parse_qual_metrics import *
import json 
import pandas as pd

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
    statsDF.to_csv(outF, sep=',') # currently = just write everything, then pick back up in plotting for subsets, etc. 
    #future: change to printing several smaller (more directed/useful) comparison files
    
def scrape_single_assembly(qualityDir,transrateDir, assembPath, d, fileD):
    ''' Extracts information from all of the relevant files in an assembly directory.
    IF ADDING ANOTHER ASSESSMENT TOOL: Add a line for the tool here; add a file parser to parse_qual_metrics.py '''
    #set up dir paths
    qualityPath = os.path.join(assembPath, qualityDir)
    transratePath = os.path.join(qualityPath, transrateDir)
    #scrape info into a dict
    infoD = {}
    for prog, fileName in fileD.items():
    ### transrate post assembly vs transrate quality --> need to figure out how to properly deal with input read stats; etc. Or just don't? That affect 'score' too..
        if prog == 'transrate': # the d from here can be diff btwn assemblies, depending on whether reads were used or not... 
    	    transrateD = transrate_parser(transratePath, fileName) 
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
	    fastqcPreTrimD = fastqc_parser(os.path.join(assembPath, fastqc_pre_trim), fileName)
	    fastqcFinalD = fastqc_parser(os.path.join(assembPath, fastqc_final), fileName)
	    infoD.update(fastqcPreTrimD)
	    infoD.update(fastqcFinalD)
    print d
    print infoD
    return infoD


def main(assembDir, qualityDir, transrateDir, outN, transratePostAssemblyDir): 
    outFile = os.path.join(assembDir, outN)  #general output filename (for both .csv and .json) 
    assembliesD = {} 
    for d in os.walk(assembDir).next()[1]: # list *only* directories
        if is_assembly_dir(assembDir,d, qualityDir, transrateDir):
            assemblyPath = os.path.join(assembDir, d)
	    assemblyName = d 
	    assembliesD[assemblyName] = scrape_single_assembly(qualityDir,transrateDir, assemblyPath, d, fileD) 
        print assembliesD
    assembliesDF = pd.DataFrame.from_dict(assembliesD)
    writeAssemblyStats_csv(assembliesDF, outFile)
    #writeAssemblyStats_json(assembliesDF, outFile)



if(__name__ == '__main__'):
    parser = argparse.ArgumentParser(description="Description: Gather stats on assemblies performed by MMT.")
    parser.add_argument('-a','--assemblies_home', help='The directory that contains the assembly directories.')
    parser.add_argument('-o','--outBase', help='output assembly stats basename. Results will print to outBase.csv and outBase.json')
    parser.add_argument('-q','--qualityDir', help='optional: alternative name for quality directory', default='quality_files')
    parser.add_argument('-t','--transrateDir', help='optional: alternative name for transrate directory', default='transrate_quality')
    parser.add_argument('--transratePostAssemblyDir', help='optional: alternative name for transrate directory', default='transrate_post_assembly')
    args = parser.parse_args()
    main(args.assemblies_home, args.qualityDir, args.transrateDir, args.outBase, args.transratePostAssemblyDir)


