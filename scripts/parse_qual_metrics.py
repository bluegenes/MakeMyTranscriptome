# Authors: Tessa Pierce, Nolan Hartwick

# script containing file parsers for results of assessment tools

#NOTE: relies on "run_busco_DBName" naming system for busco directories

import os, re, argparse
# use fadapa to parse FastQC metrics
from fadapa import Fadapa

def cegma_parser(fileList):
    '''	CEGMA PARSER
	Extracts the percent complete and the percent partial from an opened file rf.
	
	# only ever a single file!??? --> can we get rid of the text handling for mult files?
    '''
    cegmaInfo = {}
    if len(fileList) <1:
        cegmaInfo['cegma'] = 'no cegma information'
    else:
    #    for f in fileList:
    	f = fileList[0] # just ignore multiple cemga files right now ... shouldn't exist...
	fo = open(f, 'r')
	read_report = fo.read()
	cegmaInfo['%_complete'] = re.search('Complete\s+[0-9]+\s+.*?\s',read_report).group(0).split()[-1]
	cegmaInfo['%_partial'] = re.search('Partial\s+[0-9]+\s+.*?\s',read_report).group(0).split()[-1]
	fo.close()
    return cegmaInfo


def busco_parser(fileList):
    ''' BUSCO PARSER
        Extracts info from BUSCO's "short_summary_busco_DBNAME" file
        C:94%[D:40%],F:1.8%,M:3.5%,n:843o
    '''
    buscoInfo = {}
    if len(fileList) <1:
        buscoInfo['busco'] = 'no busco information'
    else:
        print fileList
        for f in fileList:
            buscoDbName = f.rsplit('_', 1)[1] # split on last occurance of '_' to get db name --> relies on "run_busco_DBName" naming system
            ret = {}
            fo = open(f, 'r')
            read_report = fo.read()
            info = re.search('C:(\d+\.?\d*)%\[D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%,n:(\d*)',read_report).groups()
            ret[buscoDbName+ '_' + '%_complete'] = info[0]
            ret[buscoDbName+ '_' +'%_duplicated'] = info[1]
            ret[buscoDbName+ '_' +'%_fragmented'] = info[2]
            ret[buscoDbName+ '_' +'%_missing'] = info[3]
            ret[buscoDbName+ '_' +'number_searched'] = info[4]
            buscoInfo.update(ret)
            fo.close()
    return buscoInfo

def detonate_parser(fileList):
    '''	DETONATE PARSER
	Extracts the score from the DETONATE .score file
    '''
    detonateD = {}
    if len(fileList) <1:
        detonateD['detonate'] = 'no detonate information'
    else:
	#for rf in fileList: --> later, need to deal with potential for multiple files...
	f = fileList[0]
	fo = open(f, 'r')
        read_score = fo.read()
	detonateD['detonate'] = re.search('Score\t[-+]?\d*\.\d+|\d+',read_score).group(0).split()[-1]
	fo.close()
    return detonateD

def transrate_parser(fileList):
    '''TRANSRATE PARSER
	Extracts relevant info from the transrate assemblies.csv file
	#this csv will be two lines long: first, headers, then values
    '''
    transrateInfo = {} 
    if len(fileList) <1:
        transrateInfo['transrate'] = 'no transrate information'
    else:
	f = fileList[0] # IF we need to deal with more files, change this to for loop (see busco parser)
        fo = open(f, 'r')
	headers = fo.readline().rstrip().split(',')
	values = fo.readline().rstrip().split(',')
	transrateInfo = dict(zip(headers, values))
	fo.close()	
    return transrateInfo


def transrate_readcount_parser(fileList):
    '''TRANSRATE PARSER 
    Extracts relevant info from the transrate assemblies.csv file
    '''
    readCount = {}
    if len(fileList) <1:
        readCount['readCount'] = 'no transrate readCount information'
    else:
        f = fileList[0] # IF we need to deal with more files, change this to for loop (see busco parser)
        fo = open(f, 'r')
        readCount['readCount'] = fo.readline().rstrip()
        fo.close()
    return readCount


def fastqc_parser(fileList):
    '''FASTQC PARSER
    Extracts info from fastqc output files --> currently relies on FADAPA parser
    '''
    fastqcD = {}
    if len(fileList) <1:
	fastqcD['fastqc'] = 'no fastqc information'
    else:
        f = fileList[0] # how to deal with paired vs unpaired files? Just deal with it in main assembly_scrape file? 
	fo = Fadapa(f) # fastqc_data.txt file..
	basicStats = fo.clean_data('Basic Statistics')
	fastqcD['Total Sequences'] = basicStats[4]
	fastqcD['Sequence length'] = basicStats[6]
	fastqcD['%GC'] = basicStats[7]
	fastqcD['Sequence Length Distribution'] = fo.clean_data('Sequence Length Distribution')
    return fastqcD






