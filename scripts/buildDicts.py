###############################################################################
###############################################################################
#
# Author	-	Andrew Walters
#
#
#
#
#
#
#
#
#
#
#
#
#
#
###############################################################################
###############################################################################

import urllib2

def createConversion(conversion_file): #conversion files from: http://www.genome.jp/linkdb/
    convDt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip().split()
        first = line[0].split(":")[1]
        second = line[1].split(":")[1]
        convDt[first] = second # generalized so any conv file will work. examples: swissprot: ko, go: ko, ko:pathway
    return convDt


def createConversionNogF(conversion_file):
    convDt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip().split("\t")
        first = line[0][3:]
        second = line[1]
        convDt[first] = second
    return convDt


def createConversionContigBlastNR(conversion_file):
    convDt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip().split("\t")
        first = line[0]
        second = line[1]
        third = line[2]
        fourth = line[3]
        convDt[first] = (second, third, fourth)
    return convDt


"""Need to account for duplicate contigs?"""
def createConversionContigClosest(conversion_file):
    convDt = {}
    conversion = open(conversion_file)
    prev = ""
    for line in conversion:
        line = line.strip().split("\t")
        if prev == line[0]:
        	continue
        first = line[0]
        second = line[12]
        third = line[14]
        convDt[first] = (second, third)
        prev = first
    return convDt


def createConversionEzPath(conversion_file):
    convDt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip().split("\t")
        if line[2] == 'original':
            first = line[0].split(":")[1]
            second = line[1].split(":")[1]
            convDt[first] = second.replace(",", ";")
    return convDt


def createConversionGOSlim(conversion_file):
    convDt = {}
    conversion = open(conversion_file)
    for line in conversion:
        if line[0:7] == "id: GO:":
            convDt[line[4:].strip()] = line[4:].strip()
    return convDt


"""idmapping.dat stuff"""
def createConversionIDMAP(conversion_file):
    convDt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.strip("\n").split("\t")
        convDt[line[0]] = line[2].replace(",", ";")
    return convDt

"""end idmapping.dat stuff"""


"""SP->GO,Entrez   3->16,28"""
def createConversoinGOandEntrez(conversion_file):
    convDt = {}
    conversion = open(conversion_file)
    for line in conversion:
        line = line.split("\t")
        temp = line[6].replace(" ", "")
        convDt[line[0]] = (temp, line[2])
    return convDt

