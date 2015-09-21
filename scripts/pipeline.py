import argparse
import os
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf
from assembler import gen_assembly_supervisor
from itertools import chain


def get_args():
	parser = argparse.ArgumentParser(description="Descritpion: main script for executing pipeline.")
	parser.add_argument('-no_log',help='Pipeline will delete log files.',action='store_true')
	parser.add_argument('-test',help='Use this flag to test the pipeline.',action='store_true')
	parser.add_argument('-cegma',help='Use this flag to run cegma as part of analysis',action='store_true')
	parser.add_argument('-force',help='Use this flag to perform a fresh run of the pipeline. All steps will be executed.',action='store_true')
	parser.add_argument('-no_rmdup',action='store_true')	
	parser.add_argument('-blast_uniref90',action='store_true')
	parser.add_argument('-blast_nr',action='store_true')
	parser.add_argument('-no_trim',action='store_true')
	parser.add_argument('-rnaspades',action='store_true')
	parser.add_argument('--fasta', help='Path to a fasta file to run the second half of pipeline on.',default=None)
	parser.add_argument('--fastq', help='A list of comma seperated paths to fastq files to run pipeline on.',default=[])
	parser.add_argument('--fastq2', help='A list of comma seperated paths to fastq files that are paired with --fastq files',default=[])
	parser.add_argument('--csv', help='A csv file specifying fastq input, fastq2 input, basenames for DE, and factors for DE.',default=None)
	parser.add_argument('--unpaired', help='A lst of comma seperated paths to unpaired fastq files to be run pipeline on.',default=[])
	parser.add_argument('--sample_info', help='Path to a file specifying information for input fastq files. This information is used during DE.')
	parser.add_argument('--model', help='An optional list of comma seperated values used to run differential expression')
	parser.add_argument('--out_dir', help='Path to the ouput location. Defaults to assemblies directory inside pipeline',default=tf.PATH_ASSEMBLIES)
	parser.add_argument('--out_name', help='The name of the output directory to be made in out_dir. If unused, name will be inherited from input file names')
	parser.add_argument('--assembly_name',help='Use to set the basename of the fasta assembly within outdir/outname.',default=tf.NAME_ASSEMBLY)
	parser.add_argument('--subsample_size',help='If greater than this number of reads (in millions) is provided, sub sample down to this number. Use 0 to signal that no subsampling should be performed', default='50')
	parser.add_argument('--cpu', help='Sets the thread cap for pipeline.',default='12')
	parser.add_argument('--busco_ref',help='Set the reference that busco will use for analysis',default='metazoa')
	parser.add_argument('--email',help='Pipeline will send emails informing you of runstate')
	return parser.parse_args()


def fastq_pair_check(fastq1,fastq2):
	count1 = 0
	with open(fsatq1) as f:
		for line in f:
			if(line[0]=='@'):
				count1+=1
	count2 = 0
	with open(fastq2) as f:
		for line in f:
			if(line[0]=='@'):
				count2+=1
	return count1==count2


def check_args(args):
	if(args.test):
		if(args.out_name==None):
			args.out_name = 'test_v2' 
		args.csv = tf.PATH_SCRIPTS+'/test_data/sample_info2.csv'
	if(args.csv!=None and args.fastq!=[]):
		raise Exception('Input must be specified using either fastq files as input or a csv file as input, not both.')
	return args


def handle_csv(args):
	f = open(args.csv)
	head = f.readline().rstrip().split(',')
	args.paired_names=[]
	args.unpaired_names = []
	if(args.model==None):
		args.model = ' '.join(head[3:])
	for line in f:
		fields = line.rstrip().split(',')
		if(fields[2]!=''):
			args.paired_names.append(fields[0])
			args.fastq.append(os.path.abspath(fields[1]))
			args.fastq2.append(os.path.abspath(fields[2]))
		else:
			args.unpaired_names.append(fields[0])
			args.unpaired.append(os.path.abspath(fields[1]))
	f.close()


def gen_sample_info(args):
	args.sample_info = tf.GEN_PATH_EXPRESSION_FILES()+'/sample_info.tsv'
	si = open(args.sample_info,'w')
	f = open(args.csv)
	for line in f:
		temp = line.split(',')
		si.write('\t'.join([temp[0]]+temp[3:]))
	si.close()
	f.close()


def setup(args):
	args.cpu = int(args.cpu)
	args.subsample_size = int(args.subsample_size)
	if(args.csv!=None):
		handle_csv(args)
	if(args.unpaired==[]):
		args.unpaired = args.fastq[len(args.fastq2):]
		args.fastq = args.fastq[:len(args.fastq2)]
	input_file = args.fasta if(args.fasta!=None) else args.fastq[0] if(args.fastq!=[]) else args.unpaired[0]
	base = os.path.split(input_file)[1]
	if('.' in base):
		base = '.'.join(base.split('.')[:-1])
	if('_' in base):
		base = '_'.join(base.split()[:-1])
	tf.PATH_ASSEMBLIES = args.out_dir
	tf.NAME_OUT_DIR = args.out_name if(args.out_name!=None) else base
	tf.NAME_ASSEMBLY = args.assembly_name if(args.assembly_name!=None) else base
	if(args.force):
		pass
	tf.build_dir_task([]).run()
	if(args.csv!=None):
		gen_sample_info(args)

def main(args):
	check_args(args)
	setup(args)
	supers = []
	if( (args.fastq!=[] or args.unpaired!=[]) and args.fasta==None):
		assembly_super = gen_assembly_supervisor(args.fastq,args.fastq2,args.unpaired,args.no_trim,args.rnaspades,args.no_rmdup,args.subsample_size)
		supers.append(assembly_super)


if(__name__=='__main__'):
	args = get_args()
	main(args)
