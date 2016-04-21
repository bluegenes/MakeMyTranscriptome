##################################################################
####     Tessa Pierce
#### https://github.com/bluegenes/
#### 5.30.2013, edited 12.5.2014 and 1.13.2015
###################################################################
"""Function: Take in a FASTA file and a list of contig names: print
two fastas: one containing the fasta entries for contigs in the list; 
the other containing the remaining contigs in the input FASTA.
"""
###################################################################

import sys, re

fastaIN = open(sys.argv[1], 'r')
name_list = open(sys.argv[2], 'r')
inListOUT= open(sys.argv[3], 'w')
notInListOUT= open(sys.argv[4], 'w')

fastaLines = [ x.strip() for x in fastaIN ]
names = [x.split('\t')[0] for x in name_list]

bp = ''
contig_name = re.compile('^(>\S*)')# match just to 1st whitespace. If want to match entire line, change (>\S*) --> (>.*) 

nameList = []
# parse and store the list of contigs you want to separate out.
for contig in names:
    if not contig.startswith('>'):
    	contig = '>' + contig
    searchName = re.match(contig_name, contig).groups()[0]
    nameList.append(searchName)

#go through the fasta file, writing contigs to the inList file or notInList file.
for line in fastaLines:
    if line.startswith('>'):
        if len(bp) > 0:
	    if name in nameList:
	        inListOUT.write(full_name + '\n' + bp + '\n')
	    else:
	    	notInListOUT.write(full_name + '\n' + bp + '\n')
	    bp = ''
	name = re.match(contig_name, line).groups()[0]
	full_name = line
    else:
	bp = bp + line

#catch the last contig
if name in nameList:
    inListOUT.write(full_name + '\n' + bp + '\n')
else:
    notInListOUT.write(full_name + '\n' + bp + '\n')

fastaIN.close()
inListOUT.close()
notInListOUT.close()
name_list.close()


