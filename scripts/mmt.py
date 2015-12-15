import argparse
import os
import sys
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf
from assembler import gen_assembly_supervisor
from annotater import gen_annotation_supervisor
from expression import gen_expression_supervisor
from quality import gen_quality_supervisor
from itertools import chain
import time

#####################____Argument_Parsers____#####################
common_parser = argparse.ArgumentParser(add_help=False)
common_parser.add_argument('-no_log', help='Pipeline will delete log files.',action='store_true')
common_parser.add_argument('-force', help='Use this flag to perform a fresh run of the pipeline. All steps will be executed regradless of what has already been performed.',action='store_true')
common_parser.add_argument('--email',help='Pipeline will send emails informing you of runstate.')
common_parser.add_argument('--out_dir', help='Path to the ouput location. Defaults to assemblies directory inside pipeline',default=tf.PATH_ASSEMBLIES)
common_parser.add_argument('-o','--out_name', help='The name of the output directory to be made in out_dir. If unused, name will be inherited from input file names')
common_parser.add_argument('-test',help='Use this flag to test the pipeline.',action='store_true')

master_parser = argparse.ArgumentParser(description=('Description: MMT is a powerful convenience tool that '
        'allows a user to run a full transcriptomics pipeline with a single command. MMT will manage the '
        'assembly of the reads, annotate the resultant transcripts, and perform a differential '
        'expression analysis of your dataset. In addition to performing the above three phases (assembly, '
        'annotation, expression) in a single command, each of these three principle phases can be run '
        'individually provided appropriate input.'))

#####################____Arguments____#####################
# support assembly input
assembly_input_parser = argparse.ArgumentParser(add_help=False)
assembly_input_parser.add_argument('-a','--assembly',help='A fasta transcriptome assembly that needs to be annotated.')
# support read input 
read_input_parser = argparse.ArgumentParser(add_help=False)
read_input_parser.add_argument('-u','--unpaired',help='A comma seperated list of unpaired fastq files.')
read_input_parser.add_argument('-1','--fastq1',help='A comma seperated list of fastq files. Each file should be paired with the same indexed file in fastq2.')
read_input_parser.add_argument('-2','--fastq2',help='A comma seperated list of fastq files. Each file should be paired with the same indexed file in fastq1.')
# support csv input
csv_input_parser = argparse.ArgumentParser(add_help=False)
csv_input_parser.add_argument('--csv', help='A CSV file specifying fastq input, basenames for DE, and factors for DE. See online documentation for details on formating the CSV file.',default=None)
#support cpu input
cpu_input_parser = argparse.ArgumentParser(add_help=False)
cpu_input_parser.add_argument('--cpu', help='Sets the process cap for execution. Default is 12. Use 0 to indicate no process cap should be used.',default=12,type=int)
#ASSEMBLER ARGS
assembler_input = argparse.ArgumentParser(add_help=False)
assembler_input.add_argument('-rnaspades',help='Use this flag to specify that assembly should be performed by rnaSPAdes rather than the default Trinity.',action='store_true')
assembler_input.add_argument('-trinity_normalization',action='store_true',help='Use this flag to use the trinity normalization option')
assembler_input.add_argument('-no_rmdup',help='Use thie flag to disable the removing duplicates portion of the pre-assembly read cleaning.',action='store_true')
assembler_input.add_argument('-no_trim',help='Use this flag to disable all trimming portions of pre-assembly read cleaning. Duplicate and low quality reads will not be removed. Subsampling will still be executed.',action='store_true')
assembler_input.add_argument('-trimmomatic',help='Use trimmomatic instead of prinseq to trime reads',action='store_true')
assembler_input.add_argument('--subsample_size',help='If greater than this number of reads (in millions) is provided, sub sample down to this number. Use 0 to signal that no subsampling should be performed. The default value is 50.', default=50,type=float)
assembler_input.add_argument('--subsample_seed',help='A seed used to initialize the random number generator used during random sampling.')
assembler_input.add_argument('--truncate',help='snip reads down to this size if longer than this size. Default is no truncations.',type=int,default=-1)
#ANNOTATION ARGS
annotation_input = argparse.ArgumentParser(add_help=False)
annotation_input.add_argument('-signalp',action='store_true',help='Use this flag to execute signalP during annotation. Only use if you have installed signalP.')
annotation_input.add_argument('-tmhmm',action='store_true',help='Use this flag to execute tmhmm during annotation. Only use if you have installed tmhmm.')
annotation_input.add_argument('-rnammer',action='store_true',help='Use this flag to execute rnammer during annotation. Only use if you have installed rnammer.')
#QUALITY ARGS
quality_input = argparse.ArgumentParser(add_help=False)
quality_input.add_argument('-cegma',help='Use this flag to run cegma as part of the annotation pipeline. Cegma is an old tool for assesing the quality of assemblies. Normal behavior of the pipeline is to use busco for assesing assemblies. Using this flag will run cegma in addition to Busco.',action='store_true')
quality_input.add_argument('--transrate_ref',help='A reference that transrate will use to evaluate the quality of your assembly.',default='')

#EXPRESSION ARGS 
expression_input = argparse.ArgumentParser(add_help=False)
expression_input.add_argument('--model', help='An optional list of comma seperated values used to run differential expression. This is particularly useful for refining Differential Expression runs as it allow you to use the same input CSV file and perform new comparisons.')

#DATABASE SELECTOR ARGS 
database_selector = argparse.ArgumentParser(add_help=False)
#need to add no-metazoa arg???? so can turn off busco-metazoa?
database_selector.add_argument('-m', '--metazoa', help = 'use metazoa BUSCO database. If used with "databases" tool, download this database.',action='store_true',default=True)
database_selector.add_argument('-e', '--eukaryota', help = 'use eukaryote BUSCO database. If used with "databases" tool, download this database.',action='store_true',default=False)
database_selector.add_argument('-v', '--vertebrata', help = 'use vertebrate BUSCO database. If used with "databases" tool, download this database.',action='store_true',default=False)
database_selector.add_argument('--arthropoda', help = 'use arthropod BUSCO database. If used with "databases" tool, download this database.',action='store_true',default=False)
database_selector.add_argument('-f', '--fungi', help = 'use fungi BUSCO database. If used with "databases" tool, download this database.',action='store_true',default=False)
database_selector.add_argument('-b', '--bacteria', help = 'use bacteria BUSCO database. If used with "databases" tool, download this database.',action='store_true',default=False)
database_selector.add_argument('-p', '--plants', help = 'use plant BUSCO database. If used with "databases" tool, download this database.',action='store_true',default=False)
database_selector.add_argument('--reinstall', help= 'download new version of all databases', default=False)

annot_database_selector = argparse.ArgumentParser(add_help=False)
annot_database_selector.add_argument('-blastplus',action='store_true',help='Use the blast+ tool suite instead of diamond to align your transcripts to the references. If used with "databases" tool, download this database.')
annot_database_selector.add_argument('-uniref90',help='Use this flag to enable the uniref-90 diamond-blast runs as part of the annotation pipeline. If used with "databases" tool, download this database.',action='store_true')
annot_database_selector.add_argument('-nr',help='Use this flag to enable the NR (non-redundant protein database) diamond-blast runs as part of the annotation pipeline. If used with "databases" tool, download this database. FYI, this takes a while.',action='store_true')

#DATABASES ARGS
#database_input = argparse.ArgumentParser(add_help=False)
#database_input.add_argument('--reinstall', help= 'download new version of all databases', default=False)
#database_input.add_argument('--cpu', help= 'cpu cap for database downloads & indexing', default=4, type=int)

###################################
subparsers = master_parser.add_subparsers(title='TOOLS', description='Tool Selector', help='Select an available module')

full_parser = subparsers.add_parser('full', parents=[common_parser, cpu_input_parser, csv_input_parser, read_input_parser, assembly_input_parser, assembler_input, annotation_input, expression_input, quality_input, database_selector, annot_database_selector], description= "Selected_tool : Full. Executing this tool will run the entire transcriptomics pipeline. This tool requires a specially formatted CSV file to describe the input. Please see the online documentation to learn how to format these files.", add_help=True)
full_parser.set_defaults(which='full')

assembly_parser = subparsers.add_parser('assembly', parents=[common_parser, cpu_input_parser, csv_input_parser, read_input_parser, assembler_input], description='Selected_tool : Assembler. Executing this tool will  clean all provided reads and assemble them.', add_help=True)
assembly_parser.set_defaults(which='assembly')

annotation_parser = subparsers.add_parser("annotation", parents=[common_parser, cpu_input_parser, assembly_input_parser, annotation_input, annot_database_selector], description='Selected_tool : Annotation. Executing this tool will run asssembly quality assesment along with a series of tools designing to provide information about the assembled transcripts.', add_help=True)
annotation_parser.set_defaults(which='annotation')

quality_parser = subparsers.add_parser("quality", parents=[common_parser, cpu_input_parser, csv_input_parser, read_input_parser, assembly_input_parser, quality_input], description='Selected_tool : Quality. Executing this tool will perform a quality assessment of an existing assembly.', add_help=True)
quality_parser.set_defaults(which='quality')

expression_parser = subparsers.add_parser("expression", parents=[common_parser, cpu_input_parser, csv_input_parser, expression_input], description='Selected_tool : Expression. Executing this tool will run a series of differential expression analyses and sumarize the output.',add_help=True)
expression_parser.set_defaults(which='expression')

database_parser = subparsers.add_parser('databases', parents=[cpu_input_parser, database_selector, annot_database_selector], description='Selected_tool : Databases. Executing this tool will check that all annotation databases are present and download if necessary. Optional: download new versions of all databases', add_help=True)
database_parser.set_defaults(which='databases')

args =  master_parser.parse_args()

#####################____Helper_Functions____#####################

def inherit_name(filePath):
    base = os.path.basename(filePath)
    base = base if('.' not in base) else '.'.join(base.split('.')[:-1])
    return base

def inherit_name_from_paired(paired_1):
    base = os.path.basename(paired_1)
    base = base if('.' not in base) else '.'.join(base.split('.')[:-1])
    base = base if('1' not in base) else '1'.join(base.split('1')[:-1])
    return base

def set_test_args(args):
    args.out_name = args.out_name if(args.out_name!=None) else 'test_v2'
    if(args.which=='full' or args.which=='assembly' or args.which=='expression' or args.which=='quality'):
        args.csv = tf.PATH_SCRIPTS+'/test_data/sample_info2.csv'
    if(args.which=='quality' or args.which=='annotation' or args.which=='expression'):
        args.assembly = tf.PATH_SCRIPTS+'/test_data/test_v2.fasta'

def set_subsample_size(args):
    args.subsample_size = 10**15 if(args.subsample_size<=0) else args.subsample_size*10**6

def gen_sample_info(args):
    args.sample_info = tf.GEN_PATH_EXPRESSION_FILES()+'/sample_info.tsv'
    si = open(args.sample_info,'w')
    f = open(args.csv)
    for line in f:
        temp = line.split(',')
        si.write('\t'.join([temp[0]]+temp[3:]))
    si.close()
    f.close()

def global_setup(args):
    if(not os.path.isdir(args.out_dir)):
        try:
            os.makedirs(args.out_dir)
        except:
            raise Exception('\n\nERROR : invalid out_dir argument. Path does not exist.')
    if(args.cpu<=0):
        args.cpu = float('inf')
    tf.PATH_ASSEMBLIES = args.out_dir
    if(args.which=='full'):
	if(args.assembly!=None):
	    base = inherit_name(args.assembly)
        elif(args.csv!=None):
            base = inherit_name(args.csv)
	elif(args.fastq1!=[]):
	    base = inherit_name_from_paired(args.fastq1[0])
	else:
	    base = inherit_name(args.unpaired[0])
    if(args.which=='assembly'):
	if(args.csv!=None):
	    base = inherit_name(args.csv)
        elif(args.fastq1!=[]):
            base = inherit_name_from_paired(args.fastq1[0])
        else:
            base = inherit_name(args.unpaired[0])
    if(args.which=='annotation' or args.which=='quality' or args.which=='expression'):
        base = inherit_name(args.assembly)
    tf.NAME_OUT_DIR = args.out_name if(args.out_name!=None) else base
    tf.NAME_ASSEMBLY = tf.NAME_OUT_DIR
    tf.build_dir_task([]).run()


#####################____Argument_Testers____#####################

def check_csv_input(args, required=True):
    if(args.csv!=None):
	if(not os.path.isfile(args.csv)):
            raise Exception('\n\nERROR : Invalid csv argument. '+args.csv+' does not exist.')
	else:
	    handle_csv(args)
    elif required:
        raise Exception('\n\nERROR : csv input is required for execution of expression module.')

def check_read_inputs(args, required=True):
    if(args.csv!=None):
        check_csv_input(args, required)
    elif(args.fastq1 !=None or args.unpaired!=None):
	handle_read_input(args)
    elif required:
        raise Exception('\n\nERROR : No input files specified. Please specify input using either csv, or --fastq1 and --fastq2 or --unpaired.')

def check_fasta_input(args):
    if(args.assembly==None):
        raise Exception('\n\nERROR : No input assembly specified. Please specify an input assembly using --assembly.')
    if(not os.path.isfile(args.assembly)):
        raise Exception('\n\nERROR : Invalid assembly argument. '+args.assembly+' does not exist.')

def fastq_pair_check(fastq1,fastq2):
    '''
    count1 = 0
    with open(fastq1) as f:
        for line in f:
            count1+=1
    count2 = 0
    with open(fastq2) as f:
        for line in f:
            count2+=1
    return count1==count2
    '''
    return True


############## Handle Input Args ##############

def handle_read_input(args):
    args.unpaired = [] if(args.unpaired==None) else args.unpaired.split(',')
    args.fastq1 = [] if(args.fastq1==None) else args.fastq1.split(',')
    args.fastq2 = [] if(args.fastq2==None) else argsfastq2.split(',')
    if(args.unpaired!=[] or args.fastq1!=[] or args.fastq2!=[] ):
        for f in chain(args.fastq1,args.fastq2,args.unpaired):
            if(not os.path.isfile(f)):
                raise Exception('\n\nERROR : Unable to find file : '+f)
        for f1,f2 in zip(args.fastq1,args.fastq2):
            if(not fastq_pair_check(f1,f2)):
                raise Exception('\n\nERROR : '+f1+' cant be paired with '+f2+'.')


def handle_csv(args):
    f = open(args.csv)
    head = f.readline().rstrip().split(',')
    args.paired_names=[]
    args.unpaired_names=[]
    args.fastq1=[]
    args.fastq2=[]
    args.unpaired=[]
    if(args.which=='expression' or args.which=='full'):
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


######################___go_functions___######################
def go_assembly(args, dep):
    return gen_assembly_supervisor(
        args.fastq1, args.fastq2, args.unpaired, dep,args.no_trim, 
    args.rnaspades, args.no_rmdup, args.subsample_size,
        args.cpu, args.subsample_seed, args.trinity_normalization,
        args.truncate, args.trimmomatic)

def go_quality(args, dep):
    busco_args = {'arthropoda': args.arthropoda, 'metazoa': args.metazoa,
                  'vertebrata': args.vertebrata, 'eukaryota': args.eukaryota,
                  'fungi': args.fungi, 'bacteria': args.bacteria,
                  'plants': args.plants}
    busco_args = [k for k in busco_args if(busco_args[k])]
    return gen_quality_supervisor(
        args.fastq1, args.fastq2, args.unpaired, dep, busco_args,
    args.cpu, args.cegma, args.transrate_ref)

def go_annotation(args, dep):
    return gen_annotation_supervisor(
        args.cpu, args.uniref90, args.nr, args.blastplus, args.signalp,
        args.tmhmm, args.rnammer, dep)

def go_expression(args, dep):
    return gen_expression_supervisor(
        args.fastq1, args.fastq2, args.paired_names, args.unpaired,
        args.unpaired_names, args.cpu, args.sample_info, args.model, dep)

def go_manage_db(args, dep, log_files=True):
    busco_args = {'arthropoda': args.arthropoda, 'metazoa': args.metazoa,
                  'vertebrata': args.vertebrata, 'eukaryota': args.eukaryota, 
                  'fungi': args.fungi, 'bacteria': args.bacteria,
                  'plants': args.plants}
    busco_args = [k for k in busco_args if(busco_args[k])]
    return tf.manage_db_task(args.reinstall, args.nr, args.uniref90, busco_args, args.blastplus, int(round(args.cpu/4)), dep, log_files)

#####################____Main_Modules____#####################
def run_full(args):
    supers = []
    deps = []
    if(args.test):
	set_test_args(args)
    if(args.assembly==None): #we run the assembly portion of pipeline
        check_read_inputs(args, True) #reads are required
        global_setup(args)
	set_subsample_size(args)
        assembly_super = go_assembly(args, [])
        supers.append(assembly_super)
	deps.append(assembly_super)
    else: #run everything else BUT assembly
	check_fasta_input(args) #fast input is required --> check that it's a proper file
	check_read_inputs(args, False) # get read input from either the csv or command-line inputs
	global_setup(args)
    manage_db = go_manage_db(args, [])
    supers.append(manage_db)
    annotation_super = go_annotation(args, deps)
    supers.append(annotation_super)
    quality_super = go_quality(args, deps)
    supers.append(quality_super)
    if(args.csv !=None): #csv is required so we have metadata
	check_csv_input(args)# since it wasn't required earlier, we need to check that it exists, is proper file. 
        gen_sample_info(args)
        expression_super = go_expression(args, deps)
        supers.append(expression_super)
    run_supers(args, supers)


def run_assembly(args):
    if(args.test):
	set_test_args(args)
    check_read_inputs(args, True)
    global_setup(args)
    set_subsample_size(args) 
    supers = []
    assembly_super = go_assembly(args, [])
    supers.append(assembly_super)
    run_supers(args, supers)


def run_quality(args):
    if(args.test):
	set_test_args(args)
    check_fasta_input(args)
    check_read_inputs(args, False)
    global_setup(args)
    supers = []
    deps = []
    if not os.path.exists(os.path.join(tf.NAME_ASSEMBLY,'.fasta')):
        cp = tf.cp_assembly_task(args.assembly, [])
        supers.append(cp)
        deps = [cp] 
    quality_super = go_quality(args, deps)
    supers.append(quality_super)
    run_supers(args, supers)


def run_annotation(args):
    if(args.test):
	set_test_args(args)
    check_fasta_input(args)
    global_setup(args)
    supers = []
    deps = []
    if not os.path.exists(os.path.join(tf.NAME_ASSEMBLY,'.fasta')):
        cp = tf.cp_assembly_task(args.assembly, [])
        supers.append(cp)
        deps = [cp]
    annotation_super = go_annotation(args, deps)
    supers.append(annotation_super)
    run_supers(args, supers)


def run_expression(args):
    if(args.test):
	set_test_args(args)
    check_fasta_input(args)
    check_csv_input(args, True)
    global_setup(args)
    gen_sample_info(args)
    supers = []
    deps = []
    if not os.path.exists(os.path.join(tf.NAME_ASSEMBLY,'.fasta')):
        cp = tf.cp_assembly_task(args.assembly, [])
        supers.append(cp)
        deps = [cp]
    expression_super = go_expression(args, deps)
    supers.append(expression_super)
    run_supers(args, supers)


def run_databases(args):
    s = go_manage_db(args, [], False)
    s.run()

def run_supers(args, supers):
    run_log = os.path.join(tf.GEN_PATH_DIR(), 'run.log')
    total = Supervisor(tasks=supers, cpu=args.cpu, force_run=args.force, log=run_log, email=args.email)
    try:
        total.run()
    except:
        build_log(args, total.task_status)
        raise
    build_log(args, total.task_status)


def build_log(args, task_status):
    master_log = os.path.join(tf.GEN_PATH_DIR(), 'master.log')
    f = open(master_log, 'a')
    f.write('##################__'+args.which+'__##################\n')
    skipped_tasks = {t: task_status[t] for t in task_status if(task_status[t]['state'] == Supervisor.STATE_SKIPPED)}
    completed_tasks = {t: task_status[t] for t in task_status if(task_status[t]['state'] == Supervisor.STATE_FINISHED)}
    failed_tasks = {t: task_status[t] for t in task_status if(task_status[t]['state'] == Supervisor.STATE_ERR)}
    removed_tasks = {t: task_status[t] for t in task_status if(task_status[t]['state'] == Supervisor.STATE_REMOVED)}
    f.write('\nThe following jobs were detected as already having been '
            'completed. They were not executed as part of this pipeline\'s '
            'execution. Output from previous runs was used instead\n')
    for t in skipped_tasks:
        f.write('\t'+t.name+'\n')
    f.write('\nThe following jobs were exectued succesfully.\n')
    for t in completed_tasks:
        f.write('\t'+t.name+' : '+completed_tasks[t]['message']+'\n')
    f.write('\nThe following jobs encountered an unexpected error during execution.\n')
    for t in failed_tasks:
        f.write('\t'+t.name+' : '+failed_tasks[t]['message']+'\n')
    f.write('\nDue to the above errors, the following jobs could not be started.\n')
    for t in removed_tasks:
        f.write('\t'+t.name+'\n')
    f.write('\n\n')
    f.close()


if(__name__ == '__main__'):
    if(args.which=='full'):
        run_full(args)
    if(args.which=='assembly'):
        run_assembly(args)
    if(args.which=='annotation'):
        run_annotation(args)
    if(args.which=='expression'):
        run_expression(args)
    if(args.which=='quality'):
        run_quality(args)
    if(args.which=='databases'):
        run_databases(args)
