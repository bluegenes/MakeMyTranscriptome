'''
filename : randomSubset.py
Author : Nolan Hartwick
Date : 4/7/15
Descritption : 	randomSubset.py is a sampling algorithm for paired fastq files.
				It is safe to import and can also be used as a script.
Call Method : 	python randomSubset.py fastqfile1 fastqfile2 linecount numbins samplesize target1 target2
				Use "python randomSubset.py -h" for help/usage description.
USAGE MESSAGE :
	usage: randomSubset.py [-h]
	                       SourceFile1 SourceFile2 Linecount Numbins SampleSize 
	                       TargetFile1 TargetFile2

	Descritption: Fastq file sampling algorithm.

	positional arguments:
	  SourceFile1  The fastq file that is paired with SourceFile2 and that we
	               desire a sample of.
	  SourceFile2  The fastq file that is paired with SourceFile1 and that we
	               desire a sample of.
	  Linecount    The number of lines in either SourceFile1 or SourceFile2. Enter 
	  			   0 if you want randomSubset.py to compute Linecount.
	  Numbins      The number of bins to divide the Fastq file into, before taking
	               subsets of each.
	  SampleSize   The desired samplesize in number of fastq entries to include in
	               sample, or in percent as a number between 0 and 1.
	  TargetFile1  An output filename. The sample from SourceFile1 writes to
	               TargetFile1.
	  TargetFile2  An output filename. The sample from SourceFile2 writes to
	               TargetFile2.

	optional arguments:
	  -h, --help   show this help message and exit
'''
import argparse
import random

'''	Generates a psuedorandom subset of paired fastq files filename1 and 
	filename2. It first seperates all entries into numbins sequential subsets. 
	A random sample is then taken from each of those subsets.
	Params:
		filename1- Path to a fastq file
		filename2- Path to the fastq file that is paired with filename1
		linecount- the linecount of filename1/2
		numbins- the number of bins to divide the files into
		samplesize- the size of the sample in number of fastq files
		filename1target- path to write sample from filename1 to
		filename2target- path to write sample from filename2 to
'''
def GenRandomizedSubset(filename1,filename2,linecount,numbins,samplesize,filename1target,filename2target):
	entries = linecount/4
	binsize = entries/numbins
	
	#divide entries into numbins equal sized subsets. i.e. entires=100, numbins = 2, ranges = [ (0,50),(50,100) ]
	ranges = [[int(round(binsize*x)),int(round(binsize*(x+1)))] for x in range(numbins)]
	if(ranges[-1][1]!=entries):
		ranges[-1][1] = entries

	#define number of samples needed from each bin
	sampleinbins = int(samplesize/numbins)
	dif = samplesize - sampleinbins*numbins
	samplevals = [sampleinbins for x in range(numbins)]
	safetoincrement = [x for x in range(numbins) if(ranges[x][1]-ranges[x][0]>sampleinbins)]
	for x in random.sample(safetoincrement,dif):
			samplevals[x]+=1

	#choose the entries to include in final subset
	ret = []
	for x in range(len(ranges)):
		t=random.sample(range(ranges[x][0],ranges[x][1]),min(samplevals[x],ranges[x][1]-ranges[x][0])  )
		t.sort()
		ret+=t

	#generate subsets
	dump(filename1,ret,filename1target)
	dump(filename2,ret,filename2target)


'''	Takes in a fastq file to read from, and a list of entries to pull and 
	write to target
'''
def dump(filename,entries,target):
	f = open(filename)
	g = open(target,'w')
	count = 0
	pos=0
	flag=False
	for line in f:
		if(count%4==0):
			if(pos>=len(entries)):
				break
			elif(count/4==entries[pos]):
				pos+=1
				flag=True
			else:
				flag=False
		if(flag):
			g.write(line)
		count+=1
	f.close()
	g.close()

'''	Checks whether filename is a valid place to read from.
'''
def rfcheck(filename):
	argparse.FileType('r')(filename)
	return filename

'''	Checks whether filname is a valid place to write to.
'''
def wfcheck(filename):
	t = argparse.FileType('w')(filename)
	return filename

'''	Computes the number of lines in filename.
'''
def lc(filename):
	f = open(filename)
	count = 0
	for line in f:
		count+=1
	return count
	f.close()

#Acts like a main method
if(__name__=='__main__'):
	parser = argparse.ArgumentParser(description="Descritpion: Fastq file sampling algorithm.")
	parser.add_argument('SourceFile1', metavar='SourceFile1', type=rfcheck, help='The fastq file that is paired with SourceFile2 and that we desire a sample of.')
	parser.add_argument('SourceFile2', metavar='SourceFile2', type=rfcheck, help='The fastq file that is paired with SourceFile1 and that we desire a sample of.')
	parser.add_argument('Linecount', metavar='Linecount', type=int, help='The number of lines in either SourceFile1 or SourceFile2. Enter 0 if you want randomSubset.py to compute Linecount.')
	parser.add_argument('Numbins', metavar='Numbins', type=int, help = 'The number of bins to divide the Fastq file into, before taking subsets of each.')	
	parser.add_argument('SampleSize', metavar='SampleSize', type=float, help='The desired samplesize in number of fastq entries to include in sample, or in percent as a number between 0 and 1.')
	parser.add_argument('TargetFile1', metavar='TargetFile1', type=wfcheck, help='An output filename. The sample from SourceFile1 writes to TargetFile1.')
	parser.add_argument('TargetFile2', metavar='TargetFile2', type=wfcheck, help='An output filename. The sample from SourceFile2 writes to TargetFile2.') 
	args = parser.parse_args()
	if(args.Linecount==0):
		args.Linecount = lc(args.SourceFile1)
	if(args.SampleSize<1):
		args.SampleSize = int(round(args.Linecount/4*args.SampleSize))
	else:
		args.SampleSize=int(args.SampleSize)
	GenRandomizedSubset(args.SourceFile1,args.SourceFile2,args.Linecount,args.Numbins,args.SampleSize,args.TargetFile1,args.TargetFile2)
