'''
'''

from tasks_v2 import Task
import os
join = os.path.join
from functools import reduce
from path_vars import *

def build_dir_task(tasks):
	'''
	'''
	trgs = [GEN_PATH_DIR(),GEN_PATH_ASSEMBLY_FILES(),GEN_PATH_ANNOTATION_FILES(),GEN_PATH_EXPRESSION_FILES(),GEN_PATH_LOGS()]
	'''
	trgs = [GEN_PATH_DIR(),GEN_PATH_ADDITIONAL_FILES(),GEN_PATH_FASTQC_OUTPUT(),
			GEN_PATH_PRINSEQ_OUTPUT(),GEN_PATH_OTHER_OUTPUT(),GEN_PATH_TRINITY_OUTPUT(),
			GEN_PATH_CEGMA_OUTPUT(),GEN_PATH_BUSCO_OUTPUT(),GEN_PATH_BLAST_OUTPUT(),
			GEN_PATH_RNAMMER_OUTPUT(),GEN_PATH_TRANSDECODER_OUTPUT()]
	'''
	cmd = ' '.join(['mkdir -p {0!s};'.format(d) for d in trgs])
	return Task(command=cmd,dependencies=tasks,targets=trgs,stdout=os.devnull,stderr=os.devnull)


def fastqc_task(fq_files,output_name,tasks):
	'''	Defines task for running fastqc. Uses GEN_PATH_DIR(), PATH_FASTQC,
		Params :
			fq_files - list of fastq files to run fastqc on
			tasks - a list of tasks that this task is dependant on.
	'''
	trgs = ['{0!s}/fastqc_{1!s}'.format(GEN_PATH_ASSEMBLY_FILES(),output_name)]
	cmd = 'mkdir {2!s}; {0!s} {1!s} --outdir {2!s}'.format(PATH_FASTQC,' '.join(fq_files),trgs[0])
	name = 'fastqc_'+output_name
	out,err = GEN_LOGS(name)
	return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def cp_assembly_task(source,tasks):
	'''	Defines task used to initialize an assembly when running on fasta 
		files. Uses GEN_PATH_DIR() and NAME_ASSEMBLY.
		Params :
			source - The path to the source fasta that should be used for analsis
			tasks - a list of tasks that this task is dependant on.
	'''
	trgs = [GEN_PATH_ASSEMBLY()]
	cmd = 'cp {0!s} {1!s}'.format(source,trgs[0])
	name = 'setting_fasta'
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name)


def prinseq_unpaired_task(input1,basename,opts,tasks):
	'''
	'''
	trgs = ['{0!s}/{1!s}_{2!s}'.format(GEN_PATH_ASSEMBLY_FILES(),basename,os.path.basename(input1))]
	cmd = ('perl {0!s} -fastq {1!s} --out_format 3 --out_good {2!s}/{3!s} --out_bad null '
			'--trim_qual_left 20 --trim_qual_right 20 --trim_qual_type min --min_len 35 '
			'--trim_tail_left 8 --trim_tail_right 8 {4!s} -log; mv {2!s}/{3!s}.fastq {5!s}'
			).format(PATH_PRINSEQ,input1,GEN_PATH_ASSEMBLY_FILES(),basename,opts,trgs[0])
	name = basename
	out,err = GEN_LOGS(name)	
	return Task(command = cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs)


def prinseq_task(input_1, input_2, basename, opts, tasks):
	'''	Defines prinseq task. Uses GEN_PATH_DIR(), PATH_PRINSEQ
		Params :
			input_1 - a list of 1/left fastq files
			input_2 - a list of 2/right fastq files
			basename - the basename for all output files
			opts - optional params for trinity task. 
			tasks = the tasks that this task is dependant on
	'''
	trgs = ['{0!s}/{1!s}_1_{2!s}'.format(GEN_PATH_ASSEMBLY_FILES(),basename,os.path.basename(input_1)),
			'{0!s}/{1!s}_2_{2!s}'.format(GEN_PATH_ASSEMBLY_FILES(),basename,os.path.basename(input_2))]
	pseudo_trgs = ['{0!s}/{1!s}_{2!s}.fastq'.format(GEN_PATH_ASSEMBLY_FILES(),basename,x) for x in range(1,3)]
	cmd = ('perl {0!s} -fastq {1!s} -fastq2 {2!s} --out_format 3 --out_good {3!s}/{4!s} '
			'--out_bad null --trim_qual_left 20 --trim_qual_right 20 --trim_qual_type min '
			'--min_len 55 --trim_tail_left 8 --trim_tail_right 8 {5!s} -log; mv {6!s} {7!s};'
			' mv {8!s} {9!s};').format( PATH_PRINSEQ, input_1, input_2, GEN_PATH_ASSEMBLY_FILES(), 
			basename, opts,pseudo_trgs[0],trgs[0],pseudo_trgs[1],trgs[1])
	name = basename
	out,err = GEN_LOGS(name)
	return Task(command = cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs)


def remove_dups_task(left,right,out_base,tasks):
	'''	Definies rmdup_fastq_paired tasks. Uses PATH_SCRIPTS, GEN_PATH_DIR()
		Params : 
			left : a set of left/1 files to have duplicates removed from
			right : a set of right/2 files to have duplicates removed from
			out_base : Basename for output files
			tasks : A set of tasks that this task is dependant on.
	'''
	trgs = ['{0!s}/{1!s}_{2!s}.fastq'.format(GEN_PATH_ASSEMBLY_FILES(),out_base,x) for x in range(1,3)]
	cmd = ('python {0!s}/rmdup_fastq_paired.py --left {1!s} --right {2!s} '
			'--left_target {3!s} --right_target {4!s}').format(PATH_SCRIPTS,
			','.join(left),','.join(right),trgs[0],trgs[1])
	name = 'remove_dups'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,name=name,stdout=out,stderr=err,targets=trgs)


def cat_task(left,right,basename,tasks):
	'''
	'''
	trgs = ['{0!s}/{1!s}_1.fastq'.format(GEN_PATH_ASSEMBLY_FILES(),basename),
			'{0!s}/{1!s}_2.fastq'.format(GEN_PATH_ASSEMBLY_FILES(),basename)]
	cmd = 'cat {0!s} > {1!s}; cat {2!s} > {3!s}'.format(' '.join(left),trgs[0],' '.join(right),trgs[1])
	name = 'cat_basename'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,name=name,stdout=out,stderr=err,targets=trgs)


def subset_task(infiles,out_base,num,tasks):
	'''	Defines subset task. Uses GEN_PATH_DIR(), PATH_SCRIPTS.
		Params :
			infiles - a pair of fastq files to read from
			out_base - basename of output files.
			num - the number of reads to keep.
			tasks - a list of tasks that this is dependant on.
	'''
	trgs = ['{0!s}/{1!s}_{2!s}.fastq'.format(GEN_PATH_ASSEMBLY_FILES(),out_base,x) for x in range(1,3)]
	cmd = 'python {0!s}/random_subset.py {1!s} 0 100 {2!s} {3!s} {4!s}'.format(
			PATH_SCRIPTS,' '.join(infiles),num,trgs[0],trgs[1])
	name = 'subset_reads'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,name=name,stdout=out,stderr=err,targets=trgs)


def trinity_task(fastq,fastq2,unpaired,cpu_cap, tasks):
	'''	Defines the trinity task. Uses GEN_PATH_DIR(), PATH_TRINITY, NAME_ASSEMBLY
		Params :	
			left - a 1/left fastq files
			right - a 2/right fastq files
			cpu_cap - number of threads used by trinity
			tasks - a list of tasks that this task is dependant on
	'''
	input_str = ''
	if(unpaired!=[] and fastq==[]):
		input_str+='--single '+','.join(unpaired)
	if(fastq!=[]):
		input_str+='--left '+','.join(fastq+unpaired)
		input_str+=' --right '+','.join(fastq2)
	trgs = [GEN_PATH_ASSEMBLY()]
	cmd = ('ulimit -s unlimited; ulimit -a; {0!s} --seqType fq {1!s} --CPU {3!s}  '
			'--JM 120G --bflyCPU {3!s} --bflyHeapSpaceMax 120G --bfly_opts "-V 10 '
			'--stderr" --output {4!s}/trinity; cp {4!s}/trinity/Trinity.fasta {5!s}; '
			).format( PATH_TRINITY, input_str, 'dummy', cpu_cap, GEN_PATH_ASSEMBLY_FILES(), trgs[0])
	name = 'trinity_assembly'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,cpu=cpu_cap,stdout=out,stderr=err)


def cegma_task(cpu_cap,tasks):
	'''	Defines the cegma task. Uses PATH_DIR, PATH_CEGMA, NAME_ASSEMBLY.
		Params :
			cpu_cap - number of threads to be used by cegma
			tasks - a list of tasks that this task is dependant on (trinity_task)
	'''
	trgs = ['{0!s}/{1!s}.completeness_report'.format(GEN_PATH_ANNOTATION_FILES(),NAME_ASSEMBLY)]
	cmd = '{0!s} -g {1!s} -v -o {3!s}/{2!s} -T {4!s}'.format(PATH_CEGMA,
			GEN_PATH_ASSEMBLY(),NAME_ASSEMBLY,GEN_PATH_ANNOTATION_FILES(),cpu_cap)
	name = 'cegma'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,cpu=cpu_cap,stdout=out,stderr=err)


def busco_task(reference_name,cpu_cap,tasks):
	'''	Defines the busco task. Uses PATH_DIR, PATH_BUSCO, PATH_BUSCO_REFERENCE
		Params :
			reference_name - Name of the reference file to be used by busco
			cpu_cap - the cpu limit to be gicen to busco.
			tasks - a list of tasks that this task is dependant on.
	'''
	trgs = ['{0!s}/run_busco_{1!s}'.format(GEN_PATH_ANNOTATION_FILES(),reference_name)]
	cmd = ('cd {0!s}; /matta1/biotools/anaconda/envs/py3k/bin/python {1!s} '
			'-o busco_{2!s} -in {3!s} -l {4!s}/{2!s} -m trans -f -c {5!s}'
			).format(GEN_PATH_ANNOTATION_FILES(),PATH_BUSCO,reference_name,GEN_PATH_ASSEMBLY(),
			PATH_BUSCO_REFERENCE,cpu_cap)
	name = 'busco_'+reference_name
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,cpu=cpu_cap,stdout=out,stderr=err)


def assembly_stats_task(tasks):
	''' Defines assembly_stats task. Uses PATH_DIR, PATH_SCRIPTS, NAME_ASSEMBLY.
		Params :
			tasks - a list of tasks that this task is dependant on (trinity_task)
	'''
	trgs = ['{0!s}/{1!s}.stats'.format(GEN_PATH_ANNOTATION_FILES(),NAME_ASSEMBLY)]
	cmd = 'python {0!s}/assembly_stats.py {1!s}/{2!s}.fasta > {3!s}'.format(PATH_SCRIPTS,GEN_PATH_DIR(),NAME_ASSEMBLY,trgs[0])
	name = 'assembly_stats'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def gene_trans_map_task(tasks):
	'''	Defines gene_trans_map task. Uses NAME_ASSEMBLY, PATH_DIR, PATH_GENE_TRANS_MAP.
		Params :
			tasks - a list of tasks that this task is dependant on (trinity_task) 
	'''
	trgs = ['{0!s}/{1!s}.gene_trans_map'.format(GEN_PATH_ANNOTATION_FILES(),NAME_ASSEMBLY)]
	cmd = '{0!s} {1!s}/{2!s}.fasta > {3!s}'.format(PATH_GENE_TRANS_MAP,GEN_PATH_DIR(),NAME_ASSEMBLY,trgs[0])
	name = 'gene_trans_map'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def blastx_task(path_db,cpu_cap,tasks):
	''' Defines blastx task. Uses PATH_DIR, NAME_ASSEMBLY, PATH_BLASTX.
		Params :
			path_db - a path to a blastx database
			output_name - name of output
			cpu_cap - number of threads used by task
			tasks - a list of tasks that this task is dependant on (trinity_task)
	'''
	db_name = os.path.basename(path_db).split('.')[0]
	trgs = ["{0!s}/{1!s}_{2!s}.blastx".format(GEN_PATH_ANNOTATION_FILES(), NAME_ASSEMBLY,db_name)]
	cmd = ('{0!s} -query {1!s}/{2!s}.fasta -db {3!s} -num_threads {4!s} -max_target_seqs 1 '
			'-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart '
			'send evalue bitscore stitle qcovs slen" -evalue 0.0001 > {5!s}'
			).format( PATH_BLASTX, GEN_PATH_DIR(), NAME_ASSEMBLY, path_db, cpu_cap, trgs[0])
	name = 'blastx_{0!s}'.format(db_name)
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,cpu=cpu_cap,stdout=out,stderr=err)

	
def rnammer_task(tasks):
	''' defines rnammer task. Uses NAME_ASSEMBLY, PATH_DIR, PATH_RNAMMER_PL, PATH_RNAMMER.
		Params :
			tasks - a list of tasks that this task is dependant on (trinity_task)
	'''
	trgs = ['{0!s}/{1!s}.fasta.rnammer.gff'.format(GEN_PATH_ANNOTATION_FILES(),NAME_ASSEMBLY)]
	cmd = ("cd {0!s}; {1!s} --transcriptome {2!s}/{3!s}.fasta --path_to_rnammer {4!s} "
			"--org_type euk; cd -").format(GEN_PATH_ANNOTATION_FILES(),PATH_RNAMMER_PL,
			GEN_PATH_DIR(),NAME_ASSEMBLY,PATH_RNAMMER)
	name = 'rnammer'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def predict_orfs_task(cpu_cap, tasks):
	'''	defines predict_orfs task. Uses NAME_ASSEMBLY, PATH_DIR, PATH_TRANSDECODER.
		Params : 
			cpu_cap - the number of threads to be used by task
			tasks - a list of tasks that this task is dependant on (trinity_task)
		blasp, pfamm, tmhmm, signalp are children

		*trasndecoder.* are targets
	'''
	trgs = ['{0!s}/{1!s}.fasta.transdecoder.pep'.format(GEN_PATH_ANNOTATION_FILES(),NAME_ASSEMBLY)]
	cmd = ("cd {0!s}; {1!s} -t {2!s} --workdir {0!s} --CPU {3!s} 2>&1 | "
			"tee {0!s}/transDecoder_no_pfam_log").format(GEN_PATH_ANNOTATION_FILES(),
			PATH_TRANSDECODER,GEN_PATH_ASSEMBLY(),cpu_cap)
	name = 'predict_orfs'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def signalp_task(pep_path, tasks):
	'''	Defines signalp task. Uses PATH_DIR, PATH_SIGNALP.
		Params : 
			pep_path - path to the pep file to run signalp on
			output_name - name of output file
			tasks - a list of tasks that this task is dependant on.
	'''
	output_name = 'signalp'
	trgs = ['{0!s}/{1!s}'.format(GEN_PATH_ANNOTATION_FILES(),output_name)]
	cmd = '{2!s} -f short -n {1!s} {0!s} > {1!s}.log'.format(pep_path,trgs[-1],PATH_SIGNALP)
	name = 'signalp'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def blastp_task(pep_path,path_db,cpu_cap,tasks):
	'''	Defines a task for running blastp. Uses PATH_DIR, PATH_BLASTP.
		Params : 
			pep_path - path to the pep file to run through blastp
			path_db - path to blastp databse
			cpu_cap - number of threads used by task
			tasks - a list of tasks that this task is dependant on (predict_orfs_task)
	'''
	db_name = os.path.basename(path_db).split('.')[0]
	trgs = ["{0!s}/{1!s}_{2!s}.blastn".format(GEN_PATH_ANNOTATION_FILES(), NAME_ASSEMBLY,db_name)]
	cmd = ('{0!s} -query {1!s} -db {2!s} -num_threads {3!s} -max_target_seqs 1 '
			'-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send '
			'evalue bitscore stitle qcovs slen" -evalue 0.001 > {4!s}'
			).format(PATH_BLASTP,pep_path,path_db,cpu_cap,trgs[0])
	name = 'blastp_{0!s}'.format(db_name)
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def tmhmm_task(pep_path, tasks):
	'''	Defines tmhmm task, Uses PATH_DIR, PATH_TMHMM.
		Params :
			pep_path - path to the pep file to run tmhmm on
			output_name - name of output file
			tasks - a list of tasks that this task is dependant on.
	'''
	trgs = ['{0!s}/{1!s}.tmhmm'.format(GEN_PATH_ANNOTATION_FILES(),NAME_ASSEMBLY)]
	cmd = 'cd {0!s}; {1!s} --short < {2!s} > {3!s}'.format(GEN_PATH_ANNOTATION_FILES(),PATH_TMHMM,pep_path,trgs[0])
	name = 'tmhmm'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def pfam_task(pep_path, cpu_cap, tasks):
	'''	Defines pfam task. Uses PATH_DIR, PATH_PFAM,PATH_PFAM_DB.
		Params : 
			pep_path - path to the pep file to run pfam on
			output_name - name of output file
			tasks - a list of tasks that this task is dependant on.
	'''
	trgs = ['{0!s}/{1!s}.pfam'.format(GEN_PATH_ANNOTATION_FILES(),NAME_ASSEMBLY)]
	cmd = '{4!s} --cpu {0!s} --domtblout {1!s} {2!s} {3!s} > {1!s}.log'.format(cpu_cap,
			trgs[0],PATH_PFAM_DB,pep_path,PATH_PFAM)
	name = 'pfam'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def annot_table_task(trans_map,blastx_sp,blastx_ur90,rnammer_results,pep_path,blastn_sp,blastn_ur90,pfam_results,signalp_results,tmhmm_results,tasks):
	'''	Defines annotTableMain task. Uses PATH_SCRIPTS, PATH_DIR, NAME_ASSEMBLY, 
		params:
			trans_map - Path to the output from gene_trans_map_task
			blastx_sp - Path to the the output from blastx_task when run against swiss prot
			blastx_ur90 - Path to the output from blastx_task when run against Uniref-90
			rnammer_results - Path to the output from rnammer_task
			pep_path - Path to the output from predict_orfs_task
			blastn_sp - Path to the output from blastn_task when run against swiss prot
			blastn_ur90 - path to the output from blastn_task when run against Uniref-90
			pfam_results - path to the output from pfam_task
			signalp_results - path to the output from signalp_task
			tmhmm_results - path to the output from tmhmm
			tasks - a list of tasks that this task is dependant on
	'''
	suffixes = ['annotation.txt','annotation_by_gene.txt']
	trgs = ['{0!s}/{1!s}_{2!s}'.format(GEN_PATH_DIR(),NAME_ASSEMBLY,sufx) for sufx in suffixes]
	cmd = ('python {0!s}/annotTableMain.py --fasta {1!s} --geneTransMap {2!s} ' 
			'--spX {3!s} --ur90X {4!s} --rnammer {5!s} --transdecoder {6!s} --spP {7!s} '
			' --ur90P {8!s} --pfam {9!s} --signalP {10!s} --tmhmm {11!s} --ko2path {12!s}/orthology_pathway.list '
			'--sp2enzyme {12!s}/swiss_enzyme.list --enzyme2path {12!s}/enzyme_pathway.list '
			'--pfam2enzyme {12!s}/pfam_enzyme.list --go2path {12!s}/go_pathway.txt '
			'--nog2function {12!s}/allKOG_functional_info.txt --contig2closest {12!s}/contig2closest '
			'--go2slim {12!s}/goslim_generic.obo --contig2blastnr {12!s}/contig2blastnr.txt '
			'--sp2ko {12!s}/idmapping.KO --sp2nog {12!s}/idmapping.eggNOG --sp2ortho '
			'{12!s}/idmapping.orthodb --sp2bioc {12!s}/idmapping.biocyc --sp2goentrez '
			'{12!s}/idmapping_selected.tab --outfile {13!s}/{14!s}').format(PATH_SCRIPTS,GEN_PATH_ASSEMBLY(),
			trans_map,blastx_sp,blastx_ur90,rnammer_results,pep_path,blastn_sp,blastn_ur90,pfam_results,
			signalp_results,tmhmm_results,PATH_DATABASES,GEN_PATH_DIR(),NAME_ASSEMBLY)
	name = 'build_annotation_table'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def keg_task(tasks):
	'''	Defines the keg task. Uses PATH_SCRIPTS, 
		Params : 
	'''
	trgs = ['{0!s}/ko01100.pdf'.format(GEN_PATH_ANNOTATION_FILES()),
			'{0!s}/ko01100_KO.txt'.format(GEN_PATH_ANNOTATION_FILES())]
	cmd = ('python {0!s}/color_pathways2.py --path ko01100 --transcriptomeKO '
			'{1!s}/uniq_ko_annots.txt --output {2!s}').format(PATH_SCRIPTS,PATH_DATABASES,GEN_PATH_ANNOTATION_FILES())
	name='draw_keg_maps'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def build_bowtie_task(tasks):
	'''
	'''
	trgs = ['{0!s}/{1!s}.1.bt2'.format(GEN_PATH_EXPRESSION_FILES(),NAME_ASSEMBLY)]
	cmd = '{0!s}/bowtie2-build --offrate 1 -f {1!s} {2!s}/{3!s}'.format(
			PATH_BOWTIE2,GEN_PATH_ASSEMBLY(),GEN_PATH_EXPRESSION_FILES(),NAME_ASSEMBLY)
	name = 'build_bowtie'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def bowtie2_unpaired_task(fastq,out_name,opt,cpu_cap,tasks):
	'''
	'''
	opts = ['-a -t --end-to-end', '-t --local']
	trgs = ['{0!s}/{1!s}.bam'.format(GEN_PATH_EXPRESSION_FILES(),out_name)]
	cmd = ('{0!s}/bowtie2 {1!s} -L {2!s} -N 1 --threads {3!s} -x {4!s}/{5!s} -U '
			'{6!s} | samtools view -Sb - > {7!s} ').format(PATH_BOWTIE2,
			opts[opt],22,cpu_cap,GEN_PATH_EXPRESSION_FILES(),NAME_ASSEMBLY,fastq,trgs[0])
	name = 'bowtie2_'+out_name
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def bowtie2_task(fastq1,fastq2,out_name,opt,cpu_cap,tasks):
	'''	
	'''
	opts = ['-a -t --end-to-end', '-t --local']
	trgs = ['{0!s}/{1!s}.bam'.format(GEN_PATH_EXPRESSION_FILES(),out_name)]
	cmd = ('{0!s}/bowtie2 {1!s} -L {2!s} -N 1 --maxins 800 --threads {3!s} -x {4!s}/{5!s} -1 '
			'{6!s} -2 {7!s} | samtools view -Sb - > {8!s} ').format(PATH_BOWTIE2,
			opts[opt],22,cpu_cap,GEN_PATH_EXPRESSION_FILES(),NAME_ASSEMBLY,fastq1,fastq2,trgs[0])
	name = 'bowtie2_'+out_name
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def express_task(bam_input,out_name,tasks):
	'''
	'''
	trgs = ['{0!s}/{1!s}.xprs'.format(GEN_PATH_EXPRESSION_FILES(),out_name)]
	cmd = 'mkdir {1!s}/{2!s}; {0!s} --output-dir {1!s}/{2!s}/ {3!s} {4!s}; mv {1!s}/{2!s}/results.xprs {5!s}'.format(
			PATH_EXPRESS,GEN_PATH_EXPRESSION_FILES(),out_name,GEN_PATH_ASSEMBLY(),bam_input,trgs[0])
	name = 'express_'+out_name
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def counts_to_table_task(express_files,out_name,flag,tasks):
	'''
	'''
	trgs = ['{0!s}/{1!s}.countsTable'.format(GEN_PATH_EXPRESSION_FILES(),out_name)]
	count_str = ' '.join(['--counts {0!s}'.format(f) for f in express_files])
	cmd = ('python {0!s}/counts_to_table2.py --out {1!s} --inDir {2!s} '
			'--outDir {2!s} {3!s} {4!s} --geneTransMap {5!s}/{6!s}.gene_trans_map').format(
			PATH_SCRIPTS,out_name,GEN_PATH_EXPRESSION_FILES(),flag,count_str,
			GEN_PATH_ANNOTATION_FILES(),NAME_ASSEMBLY)
	name = 'counts_to_table_'+flag
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def sam_sort_task(bam_file,out_name,tasks):
	'''
	'''
	trgs = ['{0!s}/{1!s}.bam'.format(GEN_PATH_EXPRESSION_FILES(),out_name)]
	cmd = 'samtools sort {0!s} {1!s}/{2!s}'.format(bam_file,GEN_PATH_EXPRESSION_FILES(),out_name)
	name = 'sam_sort_'+out_name
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def assembly_to_bed_task(tasks):
	'''
	'''
	trgs = ['{0!s}/{1!s}.bed'.format(GEN_PATH_EXPRESSION_FILES(),NAME_ASSEMBLY)]
	cmd = 'python {0!s}/fasta_to_bed_count_length.py {1!s} {2!s}'.format(
			PATH_SCRIPTS,GEN_PATH_ASSEMBLY(),trgs[0])
	name = 'fasta_to_bed'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def intersect_bed_task(bam_file,bed_reference,output_name,tasks):
	'''
	'''
	trgs = ['{0!s}/{1!s}.bed'.format(GEN_PATH_EXPRESSION_FILES(),output_name)]
	cmd = '{0!s} intersect -abam {1!s} -b {2!s} -wb -bed > {3!s}'.format(PATH_BEDTOOLS,bam_file,bed_reference,trgs[0])
	name = 'intersect_bed'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def deseq2_task(counts_to_table_results,sample_info,basename,model,tasks):
	'''
	'''
	trgs = ['{0!s}/deseq2_{1!s}_Condition/'.format(GEN_PATH_EXPRESSION_FILES(),basename)]
	cmd = 'Rscript {5!s}/deseq2.r --args {0!s} {1!s} {2!s} {3!s} {4!s}; cp {1!s} {2!s}/sample_info.txt;'.format(
			counts_to_table_results,sample_info,GEN_PATH_EXPRESSION_FILES(),basename,model,PATH_SCRIPTS)
	name = 'de_'+basename
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def rnaspades_task(left,right,unpaired,cpu_cap,tasks):
	'''
	'''
	virtual_target = '{0!s}/rna_spades_out_dir'.format(GEN_PATH_ASSEMBLY_FILES())
	trgs = [GEN_PATH_ASSEMBLY()]
	input_strings = []
	if(left!=None):
		input_strings.append('-1 '+left)
		input_strings.append('-2'+right)
	if(unpaired!=None):
		input_strings.append('-s '+unpaired)
	cmd = '{0!s} {1!s} --threads {2!s} -o {3!s}; cp {3!s}/contigs.fasta {4!s};'.format(
			PATH_RNASPADES,' '.join(input_strings),cpu_cap,virtual_target,trgs[0])
	name = 'rnaSPAdes_assembly'
	out,err = GEN_LOGS(name)
	return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


