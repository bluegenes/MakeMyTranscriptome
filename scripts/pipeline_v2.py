import argparse
import os
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf
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
	parser.add_argument('--subsample_size',help='If greater than this number of reads (in millions) is provided, sub sample down to this number. Use "none" to signal that no subsampling should be performed', default='50')
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


def build_assembly2(fastq1,fastq2,unpaired,args,dependency_set):
	tasks = []
	tasks.append(tf.fastqc_task(fastq1+fastq2+unpaired,'pre_trimming',[]))
	if(not args.no_trim):
		prinseq_count=0
		unpaired_prinseq_tasks=[]
		paired_prinseq_tasks = []
		prinseq_opts = '--derep 14' if(args.rmdup) else ''
		for input1 in unpaired:
			p_task = tf.prinseq_unpaired_task(input1,'prinseq_output_'+str(prinseq_count),prinseq_opts,[])
			prinseq_count+=1
			unpaired_prinseq_tasks.append(p_task)
		for input1,input2 in zip(fastq1,fastq2):
			p_task = tf.prinseq_task(input1,input2,'prinseq_output_'+str(prinseq_count),prinseq_opts,[])
			prinseq_count+=1
			paired_prinseq_tasks.append(p_task)
	tasks.extend(paired_prinseq_tasks)
	tasks.extend(unpaired_prinseq_tasks)
	unpaired_prinseq_targets = [t.targets[0] for t in unpaired_prinseq_tasks]
	paired_prinseq_targets = reduce(list.__add__,[t.targets for t in paired_prinseq_tasks])
	tasks.append(tf.fastqc_task(unpaired_prinseq_targets+paired_prinseq_targets,
								'post_trimming',unpaired_prinseq_tasks+paired_prinseq_tasks))
	trinity_input_left = []
	trinity_input_right = []
	trinity_dep = [t for t in unpaired_prinseq_tasks]
	if(fastq!=[]):
		prinseq_left = [t.targets[0] for t in paired_prinseq_tasks]
		prinseq_right = [t.targets[1] for t in paired_prinseq_tasks]
		if(rmdup):
			subset_dep = tf.remove_dups_task(prinseq_left,prinseq_right,'unique',paired_prinseq_tasks)
		else: 
			subset_dep = tf.cat_task(prinseq_left,prinseq_right,'joined',paired_prinseq_tasks)
		tasks.append(subset_dep)
		subset = tf.subset_task(subset_dep.targets,'trinity_input',subsample_size*10**6,[subset_dep])
		tasks.append(subset)
		trinity_dep.append(subset)
		trinity_input_left.append(subset.targets[0])
		trinity_input_right.append(subset.targets[1])
	tasks.append(tf.trinity_task(trinity_input_left,trinity_input_right,unpaired_prinseq_targets,cpu,trinity_dep))
	return Supervisor(tasks=tasks,cpu=cpu,email=email)


def build_assembly(args,dependency_set):
	if(args.no_trim):
		rmdup = tf.cat_task(args.fastq,args.fastq2,'all',[])
		subset = tf.subset_task(rmdup.targets,'trinity_input',args.subsample_size*10**6,[rmdup])
		late_fastqc = tf.fastqc_task(subset.targets+args.unpaired,'trinity_input',[subset])
		trinity = tf.trinity_task([subset.targets[0]],[subset.targets[1]],args.unpaired,
									args.cpu,[subset])
		all_tasks = [rmdup,subset,late_fastqc,trinity]
		return Supervisor(tasks=all_tasks,dependencies=dependency_set)
	pre_trimming_fastqc = tf.fastqc_task(args.fastq+args.fastq2+args.unpaired,'pre_trimming',[])
	prinseq_opts = '--derep 14' if(not args.no_rmdup) else ''
	paired_prinseq_tasks = [tf.prinseq_task(fq_1,fq_2,'prinseq_output'+str(i),prinseq_opts,[]) for i,(
		fq_1,fq_2) in enumerate(zip(args.fastq,args.fastq2))]
	p_u_t = lambda x,i : tf.prinseq_unpaired_task(x,'prinseq_output_'+str(i+len(paired_prinseq_tasks)),prinseq_opts,[])
	unpaired_prinseq_tasks = [p_u_t(input1,i) for i,input1 in enumerate(args.unpaired)]
	prinseq_left = [t.targets[0] for t in paired_prinseq_tasks]
	prinseq_right = [t.targets[1] for t in paired_prinseq_tasks]
	unpaired_targets = [t.targets[0] for t in unpaired_prinseq_tasks]
	post_trimming_fastqc = tf.fastqc_task(prinseq_left+prinseq_right+unpaired_targets,'post_trimming',paired_prinseq_tasks+unpaired_prinseq_tasks)
	if(args.fastq==[] or args.no_rmdup):
		trinity = tf.trinity_task(prinseq_left,prinseq_right,unpaired_targets,args.cpu,
									unpaired_prinseq_tasks+paired_prinseq_tasks)
		all_tasks = [pre_trimming_fastqc,post_trimming_fastqc,trinity]+paired_prinseq_tasks+unpaired_prinseq_tasks
	else:
		rmdup = tf.remove_dups_task(prinseq_left,prinseq_right,'unique',paired_prinseq_tasks)
		subset = tf.subset_task(rmdup.targets,'trinity_input',args.subsample_size*10**6,[rmdup])
		late_fastqc = tf.fastqc_task(subset.targets,'trinity_input',[subset])
		trinity = tf.trinity_task([subset.targets[0]],[subset.targets[1]],unpaired_targets,
									args.cpu,[subset]+unpaired_prinseq_tasks)
		all_tasks = [pre_trimming_fastqc,rmdup,subset,post_trimming_fastqc,trinity,late_fastqc] + paired_prinseq_tasks + unpaired_prinseq_tasks
	return Supervisor(tasks=all_tasks,dependencies=dependency_set)


def annotate_assembly(args,dependency_set):
	cegma = tf.cegma_task(args.cpu,[])
	busco = tf.busco_task(args.busco_ref,int(args.cpu/2),[])
	assembly_stats = tf.assembly_stats_task([])
	gene_trans_map = tf.gene_trans_map_task([])
	blastx_sprot = tf.blastx_task(tf.PATH_SWISS_PROT,int(args.cpu/2),[])
	rnammer = tf.rnammer_task([])
	predict_orfs = tf.predict_orfs_task(int(round(args.cpu/2)),[])
	pep_path = predict_orfs.targets[0]
	signalp = tf.signalp_task(pep_path,[predict_orfs])
	blastp_sprot = tf.blastp_task(pep_path,tf.PATH_SWISS_PROT,int(args.cpu/2),[predict_orfs])
	tmhmm = tf.tmhmm_task(pep_path,[predict_orfs])
	pfam = tf.pfam_task(pep_path,int(args.cpu/2),[predict_orfs])
	dependencies = [gene_trans_map,blastx_sprot,rnammer,predict_orfs,
					blastp_sprot,pfam,signalp,tmhmm]
	if(args.blast_uniref90):
		blastp_ur90 = tf.blastp_task(pep_path,tf.PATH_UNIREF90,int(args.cpu/2),[predict_orfs])
		blastx_ur90 = tf.blastx_task(tf.PATH_UNIREF90,int(args.cpu)/2,[])
		dependencies.append(blastx_ur90)
		dependencies.append(blastp_ur90)
		blastp_ur90_target = blastp_ur90.targets[0]
		blastx_ur90_target = blastx_ur90.targets[0]
	else:
		blastp_ur90_target = 'NONE'
		blastx_ur90_target = 'NONE'
	annot = tf.annot_table_task(gene_trans_map.targets[0],blastx_sprot.targets[0],
		blastx_ur90_target, rnammer.targets[0], predict_orfs.targets[0], 
		blastp_sprot.targets[0],blastp_ur90_target,pfam.targets[0],signalp.targets[0],
		tmhmm.targets[0],dependencies)
	keg = tf.keg_task([annot])
	all_tasks = [busco,assembly_stats,gene_trans_map,blastx_sprot,rnammer,
				predict_orfs,signalp,blastp_sprot,tmhmm,pfam,annot,keg]
	if(args.cegma):
		all_tasks.append(cegma)
	if(args.blast_uniref90):
		all_tasks.append(blastp_ur90)
		all_tasks.append(blastx_ur90)
	return Supervisor(tasks=all_tasks,dependencies=dependency_set)


def differential_expression(args,dependency_set):
	fasta_to_bed = tf.assembly_to_bed_task([])
	build_bowtie = tf.build_bowtie_task([])
	express_tasks = []
	bowtie_e_tasks = []
	bowtie_i_tasks = []
	sam_sort_tasks = []
	intersect_tasks = []
	for i in range(len(args.fastq)):
		bowtie_e = tf.bowtie2_task(args.fastq[i],args.fastq2[i],args.paired_names[i]+'_express_bt2',0,int(args.cpu/2),[build_bowtie])
		express = tf.express_task(bowtie_e.targets[0],args.paired_names[i],[bowtie_e])
		bowtie_i = tf.bowtie2_task(args.fastq[i],args.fastq2[i],args.paired_names[i]+'_intersect_bt2',1,int(args.cpu/2),[build_bowtie])
		sam_sort = tf.sam_sort_task(bowtie_i.targets[0],args.paired_names[i]+'_intersect_bt2_sorted',[bowtie_i])
		intersect_bed = tf.intersect_bed_task(sam_sort.targets[0],fasta_to_bed.targets[0],args.paired_names[i],[sam_sort,fasta_to_bed])
		bowtie_e_tasks.append(bowtie_e)
		express_tasks.append(express)
		bowtie_i_tasks.append(bowtie_i)
		sam_sort_tasks.append(sam_sort)
		intersect_tasks.append(intersect_bed)
	for i in range(len(args.unpaired)):
		bowtie_e = tf.bowtie2_unpaired_task(args.unpaired[i],args.unpaired_names[i]+'_express_bt2',0,int(args.cpu/2),[build_bowtie])
		express = tf.express_task(bowtie_e.targets[0],args.unpaired_names[i],[bowtie_e])
		bowtie_i = tf.bowtie2_unpaired_task(args.unpaired[i],args.unpaired_names[i]+'_intersect_bt2',1,int(args.cpu/2),[build_bowtie])
		sam_sort = tf.sam_sort_task(bowtie_i.targets[0],args.unpaired_names[i]+'_intersect_bt2_sorted',[bowtie_i])
		intersect_bed = tf.intersect_bed_task(sam_sort.targets[0],fasta_to_bed.targets[0],args.unpaired_names[i],[sam_sort,fasta_to_bed])
		bowtie_e_tasks.append(bowtie_e)
		express_tasks.append(express)
		bowtie_i_tasks.append(bowtie_i)
		sam_sort_tasks.append(sam_sort)
		intersect_tasks.append(intersect_bed)
	counts_to_table_express = tf.counts_to_table_task([t.targets[0] for t in express_tasks],'express_counts','--eXpress',express_tasks)
	counts_to_table_intersect = tf.counts_to_table_task([t.targets[0] for t in intersect_tasks],'bed_counts','',intersect_tasks)
	deseq2_express = tf.deseq2_task(counts_to_table_express.targets[0],args.sample_info,'express',args.model,[counts_to_table_express])
	deseq2_intersect = tf.deseq2_task(counts_to_table_intersect.targets[0],args.sample_info,'intersect',args.model,[counts_to_table_intersect])
	all_tasks = [fasta_to_bed,build_bowtie,counts_to_table_intersect,counts_to_table_express,deseq2_intersect,deseq2_express]
	all_tasks.extend(bowtie_e_tasks+bowtie_i_tasks+express_tasks+sam_sort_tasks+intersect_tasks)
	return Supervisor(tasks=all_tasks,dependencies=dependency_set)


def main(args):
	check_args(args)
	setup(args)
	supers = []
	if( (args.fastq!=[] or args.unpaired!=[]) and args.fasta==None):
		supers.append(build_assembly(args,[]))
	else:
		supers.append(tf.cp_assembly_task(args.fasta,[]))
	supers.append(annotate_assembly(args,supers[-1:]))
	if(args.csv):
		supers.append(differential_expression(args,[supers[-2]]

			))
	log_file = os.path.join(tf.GEN_PATH_DIR(),'master.log')
	total = Supervisor(tasks=supers,log=log_file,cpu=args.cpu)
	if(args.force):
		total.force_run=True
	total.run()


if(__name__=='__main__'):
	args = get_args()
	main(args)



