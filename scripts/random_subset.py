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
from itertools import chain, izip
import time


def fq_parser(fq):
	f = open(fq)
	count = 0
	ret = []
	for line in f:
		ret.append(line)
		count+=1
		if(count%4==0):
			yield ret
			ret = []


def paired_fq_parser(fq1,fq2):
	iters = [izip(fq_parser(f1),fq_parser(f2)) for f1,f2 in izip(fq1,fq2)]
	print(iters)
	for e1,e2 in chain(*iters):
		yield (e1,e2)


def GenRandomizedSubset_v2(fq1,fq2,numbins,samplesize,target1,target2,seed=None):
	entry_count = 0
	for e in paired_fq_parser(fq1,fq2):
		entry_count+=1
	if(samplesize<1 and samplesize>0):
		samplesize = samplesize*entry_count
	if(samplesize<=0):
		samplesize = entry_count*2
	p_val = float(samplesize)/entry_count
	t1 = open(target1,'w')
	t2 = open(target2,'w')
	seed_time = seed if(seed!=None) else str(time.time())
	print("Seed Used For RNG : "+str(seed_time))
	random.seed(seed_time)
	for e1,e2 in paired_fq_parser(fq1,fq2):
		if(random.random()<p_val):
			for line in e1:
				t1.write(line)
			for line in e2:
				t2.write(line)


def GenRandomizedSubset(filename1,filename2,linecount,numbins,samplesize,filename1target,filename2target):
	entries = linecount/4
	binsize = entries/numbins
	
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
	parser = argparse.ArgumentParser(description="Descritpion: Fastq file sampling script.")
	parser.add_argument('-1','--fastq1', help='a list of comma seperated fastq files paired with the files in fastq2.')
	parser.add_argument('-2','--fastq2', help='a list of comma seperated fastq files paired with the files in fastq1.')
	parser.add_argument('-n','--numbins', type=int, help = 'The number of bins to divide the Fastq entries into, before taking subsets of each.')	
	parser.add_argument('-s','--sample_size', type=float, help='The desired samplesize in number of fastq entries to include in sample, or in percent as a number between 0 and 1.')
	parser.add_argument('-t1','--target_file1', type=wfcheck, help='An output filename. The sample from SourceFile1 writes to TargetFile1.')
	parser.add_argument('-t2','--target_file2', type=wfcheck, help='An output filename. The sample from SourceFile2 writes to TargetFile2.') 
	parser.add_argument('--seed',help='The seed value to be used by the random number generator.',default=str(time.time()))
	args = parser.parse_args()
	args.fastq1 = args.fastq1.split(',')
	args.fastq2 = args.fastq2.split(',')
	GenRandomizedSubset_v2(args.fastq1,args.fastq2,args.numbins,args.sample_size,args.target_file1,args.target_file2,args.seed)

	'''
	if(args.Linecount==0):
		args.Linecount = lc(args.SourceFile1)
	if(args.SampleSize<1):
		args.SampleSize = int(round(args.Linecount/4*args.SampleSize))
	else:
		args.SampleSize=int(args.SampleSize)
	GenRandomizedSubset(args.SourceFile1,args.SourceFile2,args.Linecount,args.Numbins,args.SampleSize,args.TargetFile1,args.TargetFile2)
	'''

