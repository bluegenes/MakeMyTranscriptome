import argparse
import os
import sys
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf
from assembler import gen_assembly_supervisor
from annotater import gen_annotation_assembly
from expression import gen_expression_super
from itertools import chain


#####################____Argument_Parsers____#####################
def master_args():
	parser = argparse.ArgumentParser(description=('Description: OCT is a powerful convenience tool that '
		'allows a user to run a full transcriptomics pipeline with a single command. OCT will manage the '
		'assembly of the reads, annotate the resultant transcripts, and perform a differential '
		'expression analysis of your dataset. In addition to performing the above three phases (assembly, '
		'annotation, expression) in a single command, each of these three principle phases can be run '
		'individually provided appropriate input.') )
	parser.add_argument('tool_selector',help=('Specifies what portions of the transcriptomics '
		'pipeline should be run. Appropriate values are "full", "assembly", "annotation", or "expression". '
		'To see the arguments available for each of these tools, simply use the tool argument followed by '
		'"-h". For example. To see all arguments available to the "assembly" tool, use "assembly -h" as '
		'arguments') )
	parser.add_argument('place_holder',nargs=argparse.REMAINDER,help=argparse.SUPPRESS)
	return parser.parse_args()

def full_args():
	parser = argparse.ArgumentParser(description="Selected_tool : Full. Executing this tool will run the entire transcriptomics pipeline. This tool requires a specially formatted CSV file to describe the input. Please see the online documentaiton to learn how to format these files.")
	parser.add_argument('tool_selector',help=argparse.SUPPRESS)
	parser.add_argument('--csv', help='A CSV file specifying fastq input, basenames for DE, and factors for DE. See online documentation for details on formating the CSV file.',default=None)
	parser.add_argument('-test',help='Use this flag to test the pipeline.',action='store_true')
	parser.add_argument('-no_rmdup',help='Use thie flag to disable the removing duplicates portion of the pre-assembly read cleaning.',action='store_true')	
	parser.add_argument('-no_trim',help='Use this flag to disable all trimming portions of pre-assembly read cleaning. Duplicate and low quality reads will not be removed. Subsampling will still be executed.',action='store_true')
	parser.add_argument('--subsample_size',help='If greater than this number of reads (in millions) is provided, sub sample down to this number. Use 0 to signal that no subsampling should be performed. The deafult value is 50.', default='50')
	parser.add_argument('-rnaspades',help='Use this flag to specify that assembly should be performed by rnaSPAdes rather than the default Trinity.',action='store_true')
	parser.add_argument('-cegma',help='Use this flag to run cegma as part of the annotation pipeline. Cegma is an old tool for assesing the quality of assemblies. Normal behavior of the pipeline is to use busco for assesing assemblies. Using this flag will run cegma in addition to Busco.',action='store_true')
	parser.add_argument('-blast_uniref90',help='Use this flag to enable the uniref-90 blast runs as part of the annotation pipeline.',action='store_true')
	parser.add_argument('-blast_nr',help='Use this flag to enable the NR (non-redundant protein database) blast runs as part of the annotation pipeline. This could take a very long time to complete.',action='store_true')
	parser.add_argument('--busco_ref',help='Set the reference that busco will use for analysis',default='metazoa')
	parser.add_argument('--model', help='An optional list of space seperated values used to run differential expression. This is particularly useful for refining Differential Expression runs as it allow you to use the same input CSV file and perform new comparisons.')
	parser.add_argument('-no_log',help='Pipeline will delete log files.',action='store_true')
	parser.add_argument('-force',help='Use this flag to perform a fresh run of the pipeline. All steps will be executed regradless of what has already been performed.',action='store_true')
	parser.add_argument('--cpu', help='Sets the process cap for execution. Default is 12. Use 0 to indicate no process cap should be used.',default=12,type=int)
	parser.add_argument('--email',help='Pipeline will send emails informing you of runstate.')
	parser.add_argument('--out_dir', help='Path to the ouput location. Defaults to assemblies directory inside pipeline',default=tf.PATH_ASSEMBLIES)
	parser.add_argument('-o','--out_name', help='The name of the output directory to be made in out_dir. If unused, name will be inherited from input file names')	
	return parser.parse_args()

def assembly_args():
	parser = argparse.ArgumentParser(description='Selected_tool : Assembler. Executing this tool will clean all provided reads and assemble them.')
	parser.add_argument('tool_selector',help=argparse.SUPPRESS)
	parser.add_argument('-u','--unpaired',help='A comma seperated list of unpaired fastq files.')
	parser.add_argument('-1','--fastq1',help='A comma seperated list of fastq files. Each file should be paired with the same indexed file in fastq2.')
	parser.add_argument('-2','--fastq2',help='A comma seperated list of fastq files. Each file should be paired with the same indexed file in fastq1.')
	parser.add_argument('-rnaspades',help='Use this flag to specify that assembly should be performed by rnaSPAdes rather than the default Trinity.',action='store_true')
	parser.add_argument('-no_rmdup',help='Use thie flag to disable the removing duplicates portion of the pre-assembly read cleaning.',action='store_true')	
	parser.add_argument('-no_trim',help='Use this flag to disable all trimming portions of pre-assembly read cleaning. Duplicate and low quality reads will not be removed. Subsampling will still be executed.',action='store_true')
	parser.add_argument('--max_reads',help='If greater than this number of reads (in millions) is provided, sub sample down to this number. Use 0 to signal that no subsampling should be performed. The deafult value is 50.', default='50')	
	parser.add_argument('-no_log',help='Pipeline will delete log files.',action='store_true')
	parser.add_argument('-force',help='Use this flag to perform a fresh run of the pipeline. All steps will be executed regradless of what has already been performed.',action='store_true')
	parser.add_argument('--cpu', help='Sets the process cap for execution. Default is 12. Use 0 to indicate no process cap should be used.',default=12,type=int)
	parser.add_argument('--email',help='Pipeline will send emails informing you of runstate.')
	parser.add_argument('--out_dir', help='Path to the ouput location. Defaults to assemblies directory inside pipeline',default=tf.PATH_ASSEMBLIES)
	parser.add_argument('-o','--out_name', help='The name of the output directory to be made in out_dir. If unused, name will be inherited from input file names')	
	return parser.parse_args()

def annotation_args():
	parser = argparse.ArgumentParser(description='Selected_tool : Annotation. Executing this tool will run asssembly quality assesment along with a series of tools designing to provide information about the assembled transcripts.')
	parser.add_argument('tool_selector',help=argparse.SUPPRESS)
	parser.add_argument('-a','--assembly',help='A fasta transcriptome assembly that needs to be annotated.')
	parser.add_argument('-cegma',help='Use this flag to run cegma as part of the annotation pipeline. Cegma is an old tool for assesing the quality of assemblies. Normal behavior of the pipeline is to use busco for assesing assemblies. Using this flag will run cegma in addition to Busco.',action='store_true')
	parser.add_argument('-blast_uniref90',help='Use this flag to enable the uniref-90 blast runs as part of the annotation pipeline.',action='store_true')
	parser.add_argument('-blast_nr',help='Use this flag to enable the NR (non-redundant protein database) blast runs as part of the annotation pipeline. This could take a very long time to complete.',action='store_true')
	parser.add_argument('-no_log',help='Pipeline will delete log files.',action='store_true')
	parser.add_argument('-force',help='Use this flag to perform a fresh run of the pipeline. All steps will be executed regradless of what has already been performed.',action='store_true')
	parser.add_argument('--cpu', help='Sets the thread cap for execution. Default is 12. Use 0 to indicate no process cap should be used.',default=12,type=int)
	parser.add_argument('--email',help='Pipeline will send emails informing you of runstate.')
	parser.add_argument('--out_dir', help='Path to the ouput location. Defaults to assemblies directory inside pipeline',default=tf.PATH_ASSEMBLIES)
	parser.add_argument('-o','--out_name', help='The name of the output directory to be made in out_dir. If unused, name will be inherited from input file names')	
	return parser.parse_args()

def expression_args():
	parser = argparse.ArgumentParser(description='Selected_tool : Expression. Executing this tool will run a series of differential expression analyses and sumarize the output.')
	parser.add_argument('tool_selector',help=argparse.SUPPRESS)
	parser.add_argument('-a','--assembly',help='A fasta transcriptome assembly that needs to be annotated.')
	parser.add_argument('--csv', help='A CSV file specifying fastq input, basenames for DE, and factors for DE. See online documentation for details on formating the CSV file.',default=None)
	parser.add_argument('--model', help='An optional list of comma seperated values used to run differential expression. This is particularly useful for refining Differential Expression runs as it allow you to use the same input CSV file and perform new comparisons.')
	parser.add_argument('-no_log',help='Pipeline will delete log files.',action='store_true')
	parser.add_argument('-force',help='Use this flag to perform a fresh run of the pipeline. All steps will be executed regradless of what has already been performed.',action='store_true')
	parser.add_argument('--cpu', help='Sets the thread cap for execution. Default is 12. Use 0 to indicate no process cap should be used.',default=12,type=int)
	parser.add_argument('--email',help='Pipeline will send emails informing you of runstate.')
	parser.add_argument('--out_dir', help='Path to the ouput location. Defaults to assemblies directory inside pipeline',default=tf.PATH_ASSEMBLIES)
	parser.add_argument('-o','--out_name', help='The name of the output directory to be made in out_dir. If unused, name will be inherited from input file names')	
	return parser.parse_args()



#####################____Argument_Testers____#####################
def global_setup(args):
	if(not os.path.isdir(args.out_dir)):
		raise Exception('\n\nERROR : invalid out_dir argument. Directory does not exist.')
	if(args.cpu<=0):
		args.cpu = float('inf')
	tf.PATH_ASSEMBLIES = args.out_dir
	if(args.tool_selector=='full'):
		base = os.path.basename(args.csv)
		base = base if('.' not in base) else '.'.join(base.split('.')[:-1])
	if(args.tool_selector=='assembly'):
		if(args.unpaired!=[]):
			base = os.path.basename(args.unpaired[0])
			base = base if('.' not in base) else '.'.join(base.split('.')[:-1])
		else:
			base = os.path.basename(args.fastq1[0])
			base = base if('.' not in base) else '.'.join(base.split('.')[:-1])
			base = base if('1' not in base) else '1'.join(base.split('1')[:-1])
	if(args.tool_selector=='annotation' or args.tool_selector=='expression'):
		base = os.path.basename(args.asssembly)
		base = base if('.' not in base) else '.'.join(base.split('.')[:-1])
	tf.NAME_OUT_DIR = args.out_name if(args.out_name!=None) else base
	tf.build_dir_task([]).run()

def check_full_args(args):
	if(args.test):
		args.out_name = args.out_name if(args.out_name!=None) else 'test_v2'
		args.csv = tf.PATH_SCRIPTS+'/test_data/sample_info2.csv'
	if(args.csv==None):
		raise Exception('\n\nERROR : csv input is required for full execution of pipeline.')
	if(not os.path.isfile(args.csv)):
		raise Exception('\n\nERROR : Invalid csv argument. '+args.csv+' does not exist.')
	args.subsample_size = 10**15 if(args.subsample_size<=0) else args.subsample_size*10**6
	handle_csv(args)

def check_assembly_args(args):
	args.subsample_size = 10**15 if(args.subsample_size<=0) else args.subsample_size*10**6
	args.unpaired = [] if(args.unpaired==None) else args.unpaired.split(',')
	args.fastq1 = [] if(args.fastq1==None) else args.fastq1.split(',')
	args.fastq2 = [] if(args.fastq2==None) else args.fastq2.split(',')
	if(args.unpaired==[] and args.fastq1==[] and args.fastq2==[]):
		raise Exception('\n\nERROR : No input files specified. Please specify input using either --fastq1 and --fastq2 or --unpaired.')
	for f in chain(args.fastq1,args.fastq2,args.unpaired):
		if(not os.path.isfile(f)):
			raise Exception('\n\nERROR : Unable to find file : '+f)
	for f1,f2 in zip(args.fastq1,args.fastq2):
		if(not fastq_pair_check(f1,f2)):
			raise Exception('\n\nERROR : '+f1+' cant be paired with '+f2+'.')

def check_annotation_args(args):
	if(args.assembly==None):
		raise Exception('\n\nERROR : No input assembly specified. Please specify an input assembly using --assembly.')
	if(not os.path.isfile(args.assembly)):
		raise Exception('\n\nERROR : Invalid assembly argument. '+args.assembly+' does not exist.')

def check_expression_args(args):
	if(args.assembly==None):
		raise Exception('\n\nERROR : No input assembly specified. Please specify an input assembly using --assembly.')
	if(not os.path.isfile(args.assembly)):
		raise Exception('\n\nERROR : Invalid assembly argument. '+args.assembly+' does not exist.')
	if(args.csv==None):
		raise Exception('\n\nERROR : csv input is required for full execution of pipeline.')
	if(not os.path.isfile(args.csv)):
		raise Exception('\n\nERROR : Invalid csv argument. '+args.csv+' does not exist.')
	handle_csv(args)

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
			args.fastq1.append(os.path.abspath(fields[1]))
			args.fastq2.append(os.path.abspath(fields[2]))
		else:
			args.unpaired_names.append(fields[0])
			args.unpaired.append(os.path.abspath(fields[1]))
	f.close()
	for f in chain(args.fastq1,args.fastq2,args.unpaired):
		if(not os.path.isfile(f)):
			raise Exception('\n\nERROR : Unable to find file : '+f)
	for f1,f2 in zip(args.fastq1,args.fastq2):
		if(not fastq_pair_check(f1,f2)):
			raise Exception('\n\nERROR : '+f1+' cant be paired with '+f2+'.')

def gen_sample_info(args):
	args.sample_info = tf.GEN_PATH_EXPRESSION_FILES()+'/sample_info.tsv'
	si = open(args.sample_info,'w')
	f = open(args.csv)
	for line in f:
		temp = line.split(',')
		si.write('\t'.join([temp[0]]+temp[3:]))
	si.close()
	f.close()



#####################____Main_Modules____#####################
def run_full(args):
	check_full_args(args)
	global_setup(args)
	gen_sample_info(args)
	supers = []
	assembly_super = gen_assembly_supervisor(args.fastq,args.fastq2,args.unpaired,[],args.no_trim,
										args.rnaspades,args.no_rmdup,args.subsample_size,args.cpu)
	supers.append(assembly_super)
	annotion_super = gen_annotation_assembly(args.cpu,args.busco_ref,args.blast_uniref90,args.cegma,[assembly_super])
	supers.append(annotion_super)
	expression_super = gen_expression_super(args.fastq,args.fastq2,args.paired_names,args.unpaired,args.unpaired_names,args.cpu,args.sample_info,args.model,[assembly_super])
	supers.append(expression_super)
	run_supers(args,supers)

def run_assembly(args):
	check_assembly_args(args)
	global_setup(args)
	supers = []
	assembly_super = gen_assembly_supervisor(args.fastq,args.fastq2,args.unpaired,[],args.no_trim,
										args.rnaspades,args.no_rmdup,args.subsample_size,args.cpu)
	supers.append(assembly_super)
	run_supers(args,supers)

def run_annotation(args):
	check_annotation_args(args)
	global_setup(args)
	supers = []
	cp = tf.cp_assembly_task(args.assembly,[])
	supers.append(cp)
	annotion_super = gen_annotation_assembly(args.cpu,args.busco_ref,args.blast_uniref90,args.cegma,[cp])
	supers.append(annotion_super)
	run_supers(args,supers)

def run_expression(args):
	check_expression_args(args)
	global_setup(args)
	gen_sample_info(args)
	supers=[]
	cp = tf.cp_assembly_task(args.assembly,[])
	supers.append(cp)
	expression_super = gen_expression_super(args.fastq,args.fastq2,args.paired_names,args.unpaired,args.unpaired_names,args.cpu,args.sample_info,args.model,[cp])
	supers.append(expression_super)
	run_supers(args,supers)

def run_supers(args,supers):
	run_log = os.path.join(tf.GEN_PATH_DIR(),'master.log')
	total = Supervisor(tasks=supers,cpu=args.cpu,force_run=args.force,log=run_log)
	total.run()



if(__name__=='__main__'):
	args = master_args()
	tool=args.tool_selector.lower()
	if(tool!='a'):
		if(tool=='full'[:len(tool)]):
			args = full_args()
			args.tool_selector='full'
			run_all(args)
		if(tool=='assembly'[:len(tool)]):
			args = assembly_args()
			args.tool_selector='assembly'
			run_assembly(args)
		if(tool=='annotation'[:len(tool)]):
			args = annotation_args()
			args.tool_selector='annotation'
			run_annotation(args)
		if(tool=='expression'[:len(tool)]):
			args = expression_args()
			args.tool_selector='expression'
			run_expression(args)
	else:
		raise Exception('\n\nError : Unable to identify what tool should be executed. Valid arguments are "full", "assembly", "annotation", or "expression".')
		
	#main(args)
