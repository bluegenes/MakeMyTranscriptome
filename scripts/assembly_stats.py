'''
filename : assembly_stats.py
Author : Nolan Hartwick
Date : 4/12/15
Descritption : 	A module for computing stats for fasta files. Module is import 
				safe and can be used as a script.
Call Method :	python assembly_stats.py FastaFile
Usage Message :
	usage: assembly_stats.py [-h] FastaFile

	Descritpion: Fasta file stats generator.

	positional arguments:
	  FastaFile   The fasta file that we desire stats for.

	optional arguments:
	  -h, --help  show this help message and exit
'''


import argparse
import json

'''	Computes the mean of an iterable l.
'''
def mean(l):
	return sum(l)/len(l)

'''	Computes the median of an iterable l.
'''
def median(l):
	l = sorted(l)
	length = len(l)
	if(length%2==0):
		return (l[int(length/2-1)]+l[int(length/2)])/2
	return l[int(length/2)]

''' Computes the nx value for a list of sequence lengths.
	PARAMS:
		lens : an iterable containing lengths of sequences in fasta file.
		x : the x value for the nx compution as a decimal. n50 would use x=0.5
'''
def fasta_nx(lens,x):
	total = sum(lens)
	lens = sorted(lens,reverse=True)
	nval = total*x
	count = 0
	for entrylen in lens:
		count+=entrylen
		if(count>nval):
			return entrylen 

'''	Simple iterator over fasta file entries in rf.
	PARAMS:
		rf : an opened readable fasta file.
'''
def fa_parser(rf):
	header = rf.readline().rstrip()
	while(header!=''):
		templine = rf.readline().rstrip()
		sequence = ''
		while(templine!='' and templine[0]!='>'):
			sequence+=templine.rstrip()
			templine=rf.readline()
		yield (header,sequence)
		header = templine.rstrip()

'''	FastaStats is a function that accepts an opened fasta file f and computes 
	some statistics for the file.
	PARAMS:
		f : a readable opened fasta file
	RETURN: dictionary containing
		'median' : int (the median length of a sequence in f.)
		'mean' : int (the mean length of a sequence in f.)
		'n50' : float (the length of the sequence containing the middle base when sequences are sorted by length)
		'total' : int (the total number of bases in all sequences in f.)
		'gcpercent' : float (the percent of all sequences that is either G or C.)
		'transcripts' : int (the number of fasta entries found in f.)
'''
def FastaStats(f):
	lens = []
	line = f.readline()
	line = f.readline()
	gccount = 0
	length=0
	for header,sequence in fa_parser(f):
		#print('test')
		lens.append(len(sequence))
		gccount+=sequence.count('G')+sequence.count('C')
	ret = {}
	ret['median'] = median(lens)
	ret['mean'] = mean(lens)
	ret['n50'] = fasta_nx(lens,0.5)
	ret['total_length'] = sum(lens)
	ret['gcpercent'] = float(gccount)/ret['total_length']
	ret['num_transcripts'] = len(lens)
	return ret

#Acts like a main method
if(__name__=='__main__'):
	parser = argparse.ArgumentParser(description="Descritpion: Fasta file stats generator.")
	parser.add_argument('f', metavar='FastaFile', type=argparse.FileType('r'), help='The fasta file that we desire stats for.')
	args = parser.parse_args()
	stats = FastaStats(args.f)
	print(json.dumps(stats))
