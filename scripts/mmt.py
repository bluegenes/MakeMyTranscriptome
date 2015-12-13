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
common_parser.add_argument('--cpu', help='Sets the process cap for execution. Default is 12. Use 0 to indicate no process cap should be used.',default=12,type=int)
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
#assembly_specific
assembly_specific = argparse.ArgumentParser(add_help=False)
assembly_specific.add_argument('-u','--unpaired',help='A comma seperated list of unpaired fastq files.')
assembly_specific.add_argument('-1','--fastq1',help='A comma seperated list of fastq files. Each file should be paired with the same indexed file in fastq2.')
assembly_specific.add_argument('-2','--fastq2',help='A comma seperated list of fastq files. Each file should be paired with the same indexed file in fastq1.')
#assembly_general
assembly_general = argparse.ArgumentParser(add_help=False)
assembly_general.add_argument('-rnaspades',help='Use this flag to specify that assembly should be performed by rnaSPAdes rather than the default Trinity.',action='store_true')
assembly_general.add_argument('-trinity_normalization',action='store_true',help='Use this flag to use the trinity normalization option')
assembly_general.add_argument('-no_rmdup',help='Use thie flag to disable the removing duplicates portion of the pre-assembly read cleaning.',action='store_true')
assembly_general.add_argument('-no_trim',help='Use this flag to disable all trimming portions of pre-assembly read cleaning. Duplicate and low quality reads will not be removed. Subsampling will still be executed.',action='store_true')
assembly_general.add_argument('-trimmomatic',help='Use trimmomatic instead of prinseq to trime reads',action='store_true')
assembly_general.add_argument('--subsample_size',help='If greater than this number of reads (in millions) is provided, sub sample down to this number. Use 0 to signal that no subsampling should be performed. The default value is 50.', default=50,type=float)
assembly_general.add_argument('--subsample_seed',help='A seed used to initialize the random number generator used during random sampling.')
assembly_general.add_argument('--truncate',help='snip reads down to this size if longer than this size. Default is no truncations.',type=int,default=-1)
#annotation_specific
annotation_specific = argparse.ArgumentParser(add_help=False)
annotation_specific.add_argument('-a','--assembly',help='A fasta transcriptome assembly that will be used for computing expression levels.', required=True)
#annotation_general
annotation_general = argparse.ArgumentParser(add_help=False)
annotation_general.add_argument('-blastplus',action='store_true',help='Use the blast+ tool suite instead of diamond to align your transcripts to the references.')
annotation_general.add_argument('-uniref90',help='Use this flag to enable the uniref-90 blast runs as part of the annotation pipeline.',action='store_true')
annotation_general.add_argument('-nr',help='Use this flag to enable the NR (non-redundant protein database) blast runs as part of the annotation pipeline. This could take a very long time to complete.',action='store_true')
annotation_general.add_argument('-signalp',action='store_true',help='Use this flag to execute signalP during annotation. Only use if you have installed signalP.')
annotation_general.add_argument('-tmhmm',action='store_true',help='Use this flag to execute tmhmm during annotation. Only use if you have installed tmhmm.')
annotation_general.add_argument('-rnammer',action='store_true',help='Use this flag to execute rnammer during annotation. Only use if you have installed rnammer.')
#quality_specific
quality_specific = argparse.ArgumentParser(add_help=False)
quality_specific.add_argument('-a','--assembly',help='A fasta transcriptome assembly that needs to be annotated.', required=True)
#quality_general
quality_general = argparse.ArgumentParser(add_help=False)
quality_general.add_argument('-cegma',help='Use this flag to run cegma as part of the annotation pipeline. Cegma is an old tool for assesing the quality of assemblies. Normal behavior of the pipeline is to use busco for assesing assemblies. Using this flag will run cegma in addition to Busco.',action='store_true')
quality_general.add_argument('--busco_ref',help='Set the reference that busco will use for analysis',default='metazoa')
quality_general.add_argument('--transrate_ref',help='A reference that transrate will use to evaluate the quality of your assembly.',default='')
#expression_general
expression_general = argparse.ArgumentParser(add_help=False)
expression_general.add_argument('--model', help='An optional list of comma seperated values used to run differential expression. This is particularly useful for refining Differential Expression runs as it allow you to use the same input CSV file and perform new comparisons.')
expression_general.add_argument('--csv', help='A CSV file specifying fastq input, basenames for DE, and factors for DE. See online documentation for details on formating the CSV file.',default=None)
#database_general
database_general = argparse.ArgumentParser(add_help=False)
database_general.add_argument('--reinstall', help= 'download new version of all databases')
database_general.add_argument('-m', '--metazoa', help = 'download metazoa BUSCO database',action='store_true',default=True)
database_general.add_argument('-e', '--eukaryota', help = 'download eukaryote BUSCO database',action='store_true',default=False)
database_general.add_argument('-v', '--vertebrata', help = 'download vertebrate BUSCO database',action='store_true',default=False)
database_general.add_argument('--arthropoda', help = 'download arthropod BUSCO database',action='store_true',default=False)
database_general.add_argument('-f', '--fungi', help = 'download fungi BUSCO database',action='store_true',default=False)
database_general.add_argument('-b', '--bacteria', help = 'download bacteria BUSCO database',action='store_true',default=False)
database_general.add_argument('-p', '--plants', help = 'download plant BUSCO database',action='store_true',default=False)
database_general.add_argument('-getUniref90',help='Download Uniref-90.',action='store_true',default=False)
database_general.add_argument('-buildBlastPlus',help='build blast+ databases.',action='store_true', default=False)


subparsers = master_parser.add_subparsers(title='TOOLS', description='Tool Selector', help='Select an available module')

full_parser = subparsers.add_parser('full', parents=[common_parser, assembly_general, annotation_general, expression_general, quality_general, database_general], description=     "Selected_tool : Full. Executing this tool will run the entire transcriptomics pipeline. This tool requires a specially formatted CSV file to describe the input. Please see the   online documentation to learn how to format these files.", add_help=True)
full_parser.set_defaults(which='full')

assembly_parser = subparsers.add_parser('assembly', parents=[common_parser,assembly_specific, assembly_general], description='Selected_tool : Assembler. Executing this tool will  clean all provided reads and assemble them.',  add_help=True)
assembly_parser.set_defaults(which='assembly')

annotation_parser = subparsers.add_parser("annotation", parents=[common_parser, annotation_specific, annotation_general], description='Selected_tool : Annotation. Executing this  tool will run asssembly quality assesment along with a series of tools designing to provide information about the assembled transcripts.', add_help=True)
annotation_parser.set_defaults(which='annotation')

quality_parser = subparsers.add_parser("quality", parents=[common_parser, quality_general], description='Selected_tool : Quality. Executing this tool will perform a quality       assessment of an existing assembly.', add_help=True)
quality_parser.set_defaults(which='quality')

expression_parser = subparsers.add_parser("expression", parents=[common_parser, expression_general], description='Selected_tool : Expression. Executing this tool will run a       series of differential expression analyses and sumarize the output.',add_help=True)
expression_parser.set_defaults(which='expression')

database_parser = subparsers.add_parser('databases', parents=[common_parser, database_general], description='Selected_tool : Annotation. Executing this tool will check that all   annotation databases are present and download if necessary. Optional: download new versions of all databases', add_help=True)
database_parser.set_defaults(which='databases')

args =  master_parser.parse_args()

#####################____Argument_Testers____#####################
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
        base = os.path.basename(args.csv)
        base = base if('.' not in base) else '.'.join(base.split('.')[:-1])
    if(args.which=='assembly'):
        if(args.unpaired!=[]):
            base = os.path.basename(args.unpaired[0])
            base = base if('.' not in base) else '.'.join(base.split('.')[:-1])
        else:
            base = os.path.basename(args.fastq1[0])
            base = base if('.' not in base) else '.'.join(base.split('.')[:-1])
            base = base if('1' not in base) else '1'.join(base.split('1')[:-1])
    if(args.which=='annotation' or args.which=='quality' or args.which=='expression'):
        base = os.path.basename(args.assembly)
        base = base if('.' not in base) else '.'.join(base.split('.')[:-1])
    tf.NAME_OUT_DIR = args.out_name if(args.out_name!=None) else base
    tf.NAME_ASSEMBLY = tf.NAME_OUT_DIR
    tf.build_dir_task([]).run()


def check_full_args(args):
    if(args.test):
        args.out_name = args.out_name if(args.out_name!=None) else 'test_v2'
        args.csv = tf.PATH_SCRIPTS+'/test_data/sample_info2.csv'
    if(args.csv==None):
        raise Exception('\n\nERROR : csv input is required for full execution of pipeline.')
    if(not os.path.isfile(args.csv)):
        raise Exception('\n\nERROR : Invalid csv argument. '+args.csv+' does not exist.')
    args.subsample_size = args.subsample_size*10**6
    handle_csv(args)


def check_assembly_args(args):
    if(args.test):
        args.out_name = args.out_name if(args.out_name!=None) else 'test_v2'
        args.csv = tf.PATH_SCRIPTS+'/test_data/sample_info2.csv'
    args.subsample_size = 10**15 if(args.subsample_size<=0) else args.subsample_size*10**6
    if(args.csv):
        handle_csv(args)
    else:
        args.unpaired = [] if(args.unpaired==None) else args.unpaired.split(',')
        args.fastq1 = [] if(args.fastq1==None) else args.fastq1.split(',')
        args.fastq2 = [] if(args.fastq2==None) else args.fastq2.split(',')
        if(args.unpaired==[] and args.fastq1==[] and args.fastq2==[] and args.csv==[]):
            raise Exception('\n\nERROR : No input files specified. Please specify input using either csv, or --fastq1 and --fastq2 or --unpaired.')
        for f in chain(args.fastq1,args.fastq2,args.unpaired):
            if(not os.path.isfile(f)):
                raise Exception('\n\nERROR : Unable to find file : '+f)
        for f1,f2 in zip(args.fastq1,args.fastq2):
            if(not fastq_pair_check(f1,f2)):
                raise Exception('\n\nERROR : '+f1+' cant be paired with '+f2+'.')

def check_quality_args(args):
    if(args.test):
        args.out_name = args.out_name if(args.out_name!=None) else 'test_v2'
        args.assembly = tf.PATH_SCRIPTS+'/test_data/test_v2.fasta'
    args.unpaired = [] if(args.unpaired==None) else args.unpaired.split(',')
    args.fastq1 = [] if(args.fastq1==None) else args.fastq1.split(',')
    args.fastq2 = [] if(args.fastq2==None) else args.fastq2.split(',')
    if(args.assembly==None):
        raise Exception('\n\nERROR : No input assembly specified. Please specify an input assembly using --assembly.')
    if(not os.path.isfile(args.assembly)):
        raise Exception('\n\nERROR : Invalid assembly argument. '+args.assembly+' does not exist.')
    for f in chain(args.fastq1,args.fastq2,args.unpaired):
        if(not os.path.isfile(f)):
            raise Exception('\n\nERROR : Unable to find file : '+f)
    for f1,f2 in zip(args.fastq1,args.fastq2):
        if(not fastq_pair_check(f1,f2)):
            raise Exception('\n\nERROR : '+f1+' cant be paired with '+f2+'.')


def check_annotation_args(args):
    if(args.test):
        args.out_name = args.out_name if(args.out_name!=None) else 'test_v2'
        args.assembly = tf.PATH_SCRIPTS+'/test_data/test_v2.fasta'
    if(args.assembly==None):
        raise Exception('\n\nERROR : No input assembly specified. Please specify an input assembly using --assembly.')
    if(not os.path.isfile(args.assembly)):
        raise Exception('\n\nERROR : Invalid assembly argument. '+args.assembly+' does not exist.')


def check_expression_args(args):
    if(args.test):
        args.out_name = args.out_name if(args.out_name!=None) else 'test_v2'
        args.csv = tf.PATH_SCRIPTS+'/test_data/sample_info2.csv'
        args.assembly = tf.PATH_SCRIPTS+'/test_data/test_v2.fasta'
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


def gen_sample_info(args):
    args.sample_info = tf.GEN_PATH_EXPRESSION_FILES()+'/sample_info.tsv'
    si = open(args.sample_info,'w')
    f = open(args.csv)
    for line in f:
        temp = line.split(',')
        si.write('\t'.join([temp[0]]+temp[3:]))
    si.close()
    f.close()


######################___go_functions___######################
def go_assembly(args, dep):
    return gen_assembly_supervisor(
        args.fastq1, args.fastq2, args.unpaired, dep,args.no_trim, 
    args.rnaspades, args.no_rmdup, args.subsample_size,
        args.cpu, args.subsample_seed, args.trinity_normalization,
        args.truncate, args.trimmomatic)

def go_quality(args, dep):
    return gen_quality_supervisor(
        args.fastq1, args.fastq2, args.unpaired, dep, args.busco_ref,
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
    return tf.manage_db_task(args.reinstall, args.getNR, args.getUniref90, busco_args, int(round(args.cpu/4)), dep, log_files)

#####################____Main_Modules____#####################
def run_full(args):
    check_full_args(args)
    global_setup(args)
    gen_sample_info(args)
    supers = []
    deps = []
    manage_db = go_manage_db(args, [])
    supers.append(manage_db)
    assembly_super = go_assembly(args, [])
    supers.append(assembly_super)
    quality_super = go_quality(args, [assembly_super])
    supers.append(quality_super)
    annotation_super = go_annotation(args, [assembly_super])
    supers.append(annotation_super)
    expression_super = go_expression(args, [assembly_super])
    supers.append(expression_super)
    run_supers(args, supers)
#try to allow assembly arg...
#if args.assembly:
#    cp = tf.cp_assembly_task(args.assembly, [])
#    supers.append(cp)
#    deps = [cp]
#    else:
#        assembly_super = go_assembly(args, [])
#        supers.append(assembly_super)
#        deps = [assembly_super]
 #   quality_super = go_quality(args, deps)
 #   supers.append(quality_super)
 #   annotation_super = go_annotation(args, deps)
 #   supers.append(annotation_super)
 #   expression_super = go_expression(args, deps)
 #   supers.append(expression_super)
 #   run_supers(args, supers)


def run_assembly(args):
    check_assembly_args(args)
    global_setup(args)
    supers = []
    assembly_super = go_assembly(args, [])
    supers.append(assembly_super)
    #quality_super = go_quality(args, [assembly_super]) # now we always run quality with assembly
    #supers.append(quality_super)
    run_supers(args, supers)


def run_quality(args):
    check_quality_args(args)
    global_setup(args)
    supers = []
    cp = tf.cp_assembly_task(args.assembly, [])
    supers.append(cp)
    quality_super = go_quality(args, [cp])
    supers.append(quality_super)
    run_supers(args, supers)


def run_annotation(args):
    check_annotation_args(args)
    global_setup(args)
    supers = []
    deps = []
    #if assembly file does not exist:
    if not os.path.exists(args.assembly):
        cp = tf.cp_assembly_task(args.assembly, [])
        supers.append(cp)
        deps = [cp]
    annotation_super = go_annotation(args, deps)
    supers.append(annotation_super)
    run_supers(args, supers)


def run_expression(args):
    check_expression_args(args)
    global_setup(args)
    gen_sample_info(args)
    supers = []
    cp = tf.cp_assembly_task(args.assembly, [])
    supers.append(cp)
    expression_super = go_expression(args, [cp])
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
