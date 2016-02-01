'''
'''

from tasks_v2 import Task
import os
from os.path import join, exists
import sys
if(sys.version[0] == '3'):
    from shutil import which
else:
    from py2_which import which_python2 as which

from external_tools import PATH_ROOT, PATH_TOOLS, TOOLS_DICT
import re

''' name variables '''
NAME_ASSEMBLY = 'myassembly'
NAME_OUT_DIR = 'mmt_test_output'

''' static path variables '''
PATH_SCRIPTS = join(PATH_ROOT, 'scripts')
PATH_DATABASES = join(PATH_ROOT, 'databases')
PATH_ASSEMBLIES = join(PATH_ROOT, 'assemblies')

''' transfer to external_tools '''
PATH_TRANSDECODER = 'TransDecoder'
PATH_TRIMMOMATIC_ADAPTERS_SINGLE = join(PATH_TOOLS,'Trimmomatic-0.35/adapters/TruSeq3-SE.fa')
PATH_TRIMMOMATIC_ADAPTERS_PAIRED = join(PATH_TOOLS,'Trimmomatic-0.35/adapters/TruSeq3-PE.fa')

''' set these up as executables that we just look for in the $path '''
PATH_BEDTOOLS = 'bedtools'
PATH_BLASTP = 'blastp'
PATH_BLASTX = 'blastx'
PATH_BOWTIE2 = ''
PATH_CEGMA = 'cegma'
PATH_EXPRESS = 'express'
PATH_KALLISTO = 'kallisto'

'''optional tools -- set up the same way as above? just don't check for them with mmt setup; only install instructions in tools''' 
PATH_RNAMMER = '/matta1/biotools/redhat/rnammer-1.2/rnammer'
PATH_RNAMMER_PL = 'RnammerTranscriptome.pl'
PATH_SIGNALP = 'signalp'
PATH_RNASPADES = 'rnaspades.py'
PATH_TMHMM = 'tmhmm'

''' static db variables '''
PATH_BUSCO_REFERENCE = join(PATH_DATABASES, 'busco')
PATH_PFAM_DATABASE = '{0!s}/pfam/Pfam-A.hmm'.format(PATH_DATABASES)
PATH_NR = os.path.join(PATH_DATABASES, 'nr', 'nr')
PATH_SWISS_PROT = os.path.join(PATH_DATABASES, 'uniprot_sprot', 'uniprot_sprot')
PATH_UNIREF90 = os.path.join(PATH_DATABASES, 'uniref90', 'uniref90')
PATH_NOG_CATEGORIES = os.path.join(PATH_DATABASES, 'nog_categories')


# Dynamic path variable functions
def GEN_PATH_DIR(): return os.path.join(PATH_ASSEMBLIES, NAME_OUT_DIR)

#def GEN_ASSEMBLY_NAME(): return NAME_ASSEMBLY

def GEN_PATH_ASSEMBLY_FILES(): return os.path.join(GEN_PATH_DIR(), 'assembly_files')

def GEN_PATH_QUALITY_FILES(): return os.path.join(GEN_PATH_DIR(), 'quality_files')

def GEN_PATH_ANNOTATION_FILES(): return os.path.join(GEN_PATH_DIR(), 'annotation_files')

def GEN_PATH_EXPRESSION_FILES(): return os.path.join(GEN_PATH_DIR(), 'expression_files')

def GEN_PATH_FILTER_FILES(): return os.path.join(GEN_PATH_DIR(), 'filtered_assemblies')

def GEN_PATH_LOGS(): return os.path.join(GEN_PATH_DIR(), 'log_files')

def GEN_PATH_ASSEMBLY(): return os.path.join(GEN_PATH_DIR(), NAME_ASSEMBLY+'.fasta')

def GEN_PATH_TRANSDECODER_DIR(): return os.path.join(GEN_PATH_ANNOTATION_FILES(), 'transdecoder')

def GEN_PATH_TRANSRATE_DIR(): return os.path.join(GEN_PATH_QUALITY_FILES(), 'transrate')

def GEN_PATH_PEP(): return os.path.join(GEN_PATH_TRANSDECODER_DIR(), NAME_ASSEMBLY+'.fasta.transdecoder.pep')

def GEN_PATH_ANNOT_TABLE(): return os.path.join(GEN_PATH_DIR(), NAME_ASSEMBLY+'annotation.txt')



def GEN_LOGS(x): return (os.path.join(GEN_PATH_LOGS(), x+'.out_log'),
                         os.path.join(GEN_PATH_LOGS(), x+'.err_log'))


def tool_path_check(name):
    temp = os.path.join(PATH_TOOLS, name)
    if(os.path.exists(temp)):
        return temp
    elif(which(name)):
        print('WARNING : MMT has not installed '+name+'. MMT will use the version found in your path variable instead.')
        return name
    else:
        #raise EnvironmentError('MMT has not installed and was unable to detect '+name)
        print('INSTALLATION ERROR : MMT has not installed '+name+' and did not find this program in your path variable. Steps requiring '+name+'will not be performed.') 

def build_dir_task(tasks):
    '''
    '''
    trgs = [GEN_PATH_DIR(), GEN_PATH_ASSEMBLY_FILES(), GEN_PATH_QUALITY_FILES(), GEN_PATH_ANNOTATION_FILES(),
            GEN_PATH_FILTER_FILES(), GEN_PATH_EXPRESSION_FILES(), GEN_PATH_LOGS()]
    cmd = ' '.join(['mkdir -p {0!s};'.format(d) for d in trgs])
    return Task(command=cmd,dependencies=tasks,targets=trgs,stdout=os.devnull,stderr=os.devnull)


def cp_assembly_task(path_assembly, source, tasks):
    '''    Defines task used to initialize an assembly when running on fasta
        files. Uses GEN_PATH_DIR() and NAME_ASSEMBLY.
        Params :
            source - The path to the source fasta that should be used for analsis
            tasks - a list of tasks that this task is dependant on.
    '''
    trgs = ['{0!s}'.format(path_assembly)]
    cmd = 'cp {0!s} {1!s}'.format(source, trgs[0]) 
    name = 'setting_fasta_' + os.path.basename(path_assembly)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name)


def fastqc_task(fq_files, output_name, tasks):
    '''    Defines task for running fastqc. Uses GEN_PATH_DIR(), PATH_FASTQC,
        Params :
            fq_files - list of fastq files to run fastqc on
            tasks - a list of tasks that this task is dependant on.
    '''
    trgs = ['{0!s}/fastqc_{1!s}'.format(GEN_PATH_ASSEMBLY_FILES(), output_name)]
    cmd = 'mkdir {2!s}; {0!s} {1!s} --outdir {2!s}'.format(
           TOOLS_DICT['fastqc'].full_exe[0],' '.join(fq_files), trgs[0])
    name = 'fastqc_'+output_name
    out, err = GEN_LOGS(name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def prinseq_unpaired_task(input1, basename, opts, tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}_{2!s}'.format(
        GEN_PATH_ASSEMBLY_FILES(), basename, os.path.basename(input1))]
    cmd = ('perl {0!s} -fastq {1!s} --out_format 3 --out_good {2!s}/{3!s} --out_bad null '
           '--trim_qual_left 20 --trim_qual_right 20 --trim_qual_type min --min_len 35 '
           '--trim_tail_left 8 --trim_tail_right 8 {4!s} -log; mv {2!s}/{3!s}.fastq {5!s}'
           ).format(tool_path_check(TOOLS_DICT['prinseq'].full_exe[0]), input1, GEN_PATH_ASSEMBLY_FILES(), 
                    basename, opts, trgs[0])
    name = basename
    out, err = GEN_LOGS(name)
    return Task(command = cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs)


def prinseq_task(input_1, input_2, basename, opts, tasks):
    '''    Defines prinseq task. Uses GEN_PATH_DIR(), PATH_PRINSEQ
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
            ' mv {8!s} {9!s};').format(tool_path_check(TOOLS_DICT['prinseq'].full_exe[0]), input_1, input_2, GEN_PATH_ASSEMBLY_FILES(), 
            basename, opts,pseudo_trgs[0],trgs[0],pseudo_trgs[1],trgs[1])
    name = basename
    out,err = GEN_LOGS(name)
    return Task(command = cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs)


def trimmomatic_unpaired_task(input1, cpu_cap, basename, tasks):
    form = lambda s, i : s.format(GEN_PATH_ASSEMBLY_FILES(), basename, os.path.basename(i))
    trgs = [form('{0!s}/{1!s}_{2!s}', input1)]
    orphans = [form('{0!s}/{1!s}_orphans_{2!s}', input1)]
    cmd = ('java -jar {0!s} SE -threads {4!s} {1!s} {2!s} {3!s} ILLUMINACLIP:'
           '{5!s}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35'
           #).format(tool_path_check(PATH_TRIMMOMATIC), input1, trgs[0], orphans[0],cpu_cap,
           ).format(tool_path_check(TOOLS_DICT['trimmomatic'].full_exe[0]), input1, trgs[0], orphans[0],cpu_cap,
           PATH_TRIMMOMATIC_ADAPTERS_SINGLE)


def trimmomatic_task(left, right, cpu_cap, basename, tasks):
    form = lambda s, i : s.format(GEN_PATH_ASSEMBLY_FILES(), basename, os.path.basename(i))
    trgs = [form('{0!s}/{1!s}_1_{2!s}', left),
            form('{0!s}/{1!s}_2_{2!s}', right)]
    orphans = [form('{0!s}/{1!s}_1s_{2!s}', left),
               form('{0!s}/{1!s}_2s_{2!s}', right)]
    cmd = ('java -jar {0!s} PE -threads {3!s} {1!s} {2!s} {5!s} {4!s} {7!s} '
           '{6!s} ILLUMINACLIP:{8!s}:2:30:10 LEADING:3 TRAILING:3 '
           'SLIDINGWINDOW:4:15 MINLEN:35').format(
           #tool_path_check(PATH_TRIMMOMATIC), left, right, cpu_cap, orphans[0], trgs[0],
           tool_path_check(TOOLS_DICT['trimmomatic'].full_exe[0]), left, right, cpu_cap, orphans[0], trgs[0],
           orphans[1], trgs[1], PATH_TRIMMOMATIC_ADAPTERS_PAIRED)
    name = basename
    out, err = GEN_LOGS(name)
    return Task(command=cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs, cpu=cpu_cap) 


def remove_dups_task(left, right, out_base, tasks):
    '''    Definies rmdup_fastq_paired tasks. Uses PATH_SCRIPTS, GEN_PATH_DIR()
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


def cat_task(left, right, basename, tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}_1.fastq'.format(GEN_PATH_ASSEMBLY_FILES(),basename),
            '{0!s}/{1!s}_2.fastq'.format(GEN_PATH_ASSEMBLY_FILES(),basename)]
    cmd = 'cat {0!s} > {1!s}; cat {2!s} > {3!s}'.format(' '.join(left),trgs[0],' '.join(right),trgs[1])
    name = 'cat_basename'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,name=name,stdout=out,stderr=err,targets=trgs)


def subset_task(fastq1, fastq2, out_base, num, seed, tasks):
    '''    Defines subset task. Uses GEN_PATH_DIR(), PATH_SCRIPTS.
        Params :
            infiles - a pair of fastq files to read from
            out_base - basename of output files.
            num - the number of reads to keep.
            tasks - a list of tasks that this is dependant on.
    '''
    trgs = ['{0!s}/{1!s}_{2!s}.fastq'.format(GEN_PATH_ASSEMBLY_FILES(),out_base,x) for x in (1,2)]
    cmd = 'python {0!s}/random_subset.py -1 {1!s} -2 {2!s} -n 100 -s {3!s} -t1 {4!s} -t2 {5!s}'.format(
            PATH_SCRIPTS,','.join(fastq1),','.join(fastq2),num,trgs[0],trgs[1])
    if(seed!=None):
        cmd+=' --seed '+seed
    name = 'subset_reads'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,name=name,stdout=out,stderr=err,targets=trgs)


def truncate_task(left, right, length, tasks):
    trgs = ['{0!s}/truncated_1.fastq'.format(GEN_PATH_ASSEMBLY_FILES(), left),
            '{0!s}/truncated_2.fastq'.format(GEN_PATH_ASSEMBLY_FILES(), right)]
    cmd = ('python {0!s}/truncate_fastq.py {1!s} --length {2!s} --target {3!s}; '
           'python {0!s}/truncate_fastq.py {4!s} --length {2!s} --target {5!s};').format(
           PATH_SCRIPTS, left, length, trgs[0], right, trgs[1])
    name = 'truncate_reads'
    out, err = GEN_LOGS(name)
    return Task(command=cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs)


def trinity_task(fastq, fastq2, unpaired, cpu_cap_trin, cpu_cap_bfly, mem_trin, mem_bfly, normalize_flag, tasks):
    '''    Defines the trinity task. Uses GEN_PATH_DIR(), PATH_TRINITY, NAME_ASSEMBLY
        Params :    
            left - a 1/left fastq files
            right - a 2/right fastq files
            cpu_cap - number of threads used by trinity
            tasks - a list of tasks that this task is dependant on
    '''
    normalize_flag = '--normalize_reads' if(normalize_flag) else ''
    input_str = ''
    if(unpaired!=[] and fastq==[]):
        input_str+='--single '+','.join(unpaired)
    if(fastq!=[]):
        input_str+='--left '+','.join(fastq+unpaired)
        input_str+=' --right '+','.join(fastq2)
    trgs = [GEN_PATH_ASSEMBLY()]
    cmd = ('{0!s} --seqType fq {1!s} --CPU {2!s} --max_memory {3!s}G --bflyCalculateCPU {4!s} '
            '--output {6!s}/trinity; cp {6!s}/trinity/Trinity.fasta {7!s};'
            #).format(tool_path_check(PATH_TRINITY), input_str, cpu_cap_trin, mem_trin, normalize_flag,
            ).format(tool_path_check(TOOLS_DICT['trinity'].full_exe[0]), input_str, cpu_cap_trin, mem_trin, normalize_flag,
                     mem_bfly, GEN_PATH_ASSEMBLY_FILES(), trgs[0])
    name = 'trinity_assembly'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,cpu=max(cpu_cap_trin,cpu_cap_bfly),stdout=out,stderr=err)


def rnaspades_task(left, right, unpaired, cpu_cap, tasks):
    '''
    '''
    virtual_target = '{0!s}/rna_spades_out_dir'.format(GEN_PATH_ASSEMBLY_FILES())
    trgs = [GEN_PATH_ASSEMBLY()]
    input_strings = []
    if(left!=[]):
        input_strings.append('-1 '+left[0])
        input_strings.append('-2 '+right[0])
    if(unpaired!=[]):
        input_strings.append('-s '+unpaired[0])
    cmd = '{0!s} {1!s} --threads {2!s} -o {3!s}; cp {3!s}/contigs.fasta {4!s};'.format(
            PATH_RNASPADES,' '.join(input_strings),cpu_cap,virtual_target,trgs[0])
    name = 'rnaSPAdes_assembly'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def cegma_task(out_dir,assembly,cpu_cap, tasks):
    '''    Defines the cegma task. Uses PATH_DIR, PATH_CEGMA, NAME_ASSEMBLY.
        Params :
            cpu_cap - number of threads to be used by cegma
            tasks - a list of tasks that this task is dependant on (trinity_task)
    '''
    assembly_name = os.path.basename(assembly).split('.fa')[0]
    trgs = ['{0!s}/{1!s}.completeness_report'.format(out_dir,assembly_name)]
    cmd = '{0!s} -g {1!s} -v -o {3!s}/{2!s} -T {4!s}'.format(PATH_CEGMA,
            assembly,assembly_name,out_dir,cpu_cap)
    name = 'cegma'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,cpu=cpu_cap,stdout=out,stderr=err)


def busco_task(assembly_path, assembly_name, out_dir,reference_name, cpu_cap, tasks):
    ''' Defines the busco task. Uses PATH_DIR, PATH_BUSCO, PATH_BUSCO_REFERENCE
        Params :
            reference_name - Name of the reference file to be used by busco
            cpu_cap - the cpu limit to be gicen to busco.
            tasks - a list of tasks that this task is dependant on.
    '''
    trgs = ['{0!s}/run_busco_{1!s}_{2!s}'.format(out_dir,assembly_name,reference_name)]
    cmd = ('cd {0!s}; /matta1/biotools/anaconda/envs/py3k/bin/python {1!s} '
            '-o busco_{3!s}_{2!s} -in {4!s} -l {5!s}/{2!s}_buscos/{2!s} -m trans -f -c {6!s}'
            ).format(out_dir,tool_path_check(TOOLS_DICT['busco_plant'].full_exe[0]),reference_name,assembly_name,assembly_path,
            PATH_BUSCO_REFERENCE,cpu_cap)
    name = 'busco_'+ reference_name + '_' + assembly_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,cpu=cpu_cap,stdout=out,stderr=err)


def transrate_dep_generator(transrate_task, lefts, rights, singles, reference, assembly_path, cpu_cap, transrate_dir, other_dependencies):

    def ret():
        for t in other_dependencies:
             try:
                if( not d.finished()):
                    return False
            except Task.ExitCodeException:
                return False
        assembly_files = sorted(os.listdir(GEN_PATH_ASSEMBLY_FILES()))
        new_lefts = [[g for g in assembly_files if(re.search(os.path.basename(f), g))][0] for f in lefts]
        new_rights = [[g for g in assembly_files if(re.search(os.path.basename(f), g))][0] for f in rights]
        new_singles = [[g for g in assembly_files if(re.search(os.path.basename(f), g))][0] for f in singles]
        if(len(new_lefts) == len(lefts) and len(new_rights) == len(rights) and len(new_singles) == len(singles)):
            lefts = ','.join(new_lefts+new_singles)
            rights = ','.join(new_rights) 
            lefts = '--left '+lefts if(len(lefts) > 0) else ''
            rights = '--right '+rights if(len(rights) > 0) else ''
            reference = '--reference ' + reference if(reference != '') else ''
            cmd = '{0!s} --assembly {1!s} {2!s} {3!s} --threads {4!s} {5!s} --output {6!s}'.format(
                   tool_path_check(TOOLS_DICT['transrate'].full_exe[0]), assembly_path, lefts,
                   rights, cpu_cap, reference, transrate_dir)
            
        else:
            print('Unable to match input files with trimmed output. Continuing transrate using input files instead.')
        return True




def transrate_task(assembly_path, assembly_name,lefts, rights, singles, out_dir, transrate_dir, cpu_cap, tasks, reference = ''): #, cpu_cap, tasks):
    trgs = ['{0!s}/assemblies.csv'.format(transrate_dir),'{0!s}/{1!s}/good.{1!s}.fasta'.format(transrate_dir,assembly_name),'{0!s}/{1!s}/{1!s}.fasta_quant.sf'.format(transrate_dir,assembly_name)]
    lefts = ','.join(lefts+singles)
    rights = ','.join(rights) 
    lefts = '--left '+lefts if(len(lefts) > 0) else ''
    rights = '--right '+rights if(len(rights) > 0) else ''
    reference = '--reference ' + reference if(reference != '') else ''
    cmd = '{0!s} --assembly {1!s} {2!s} {3!s} --threads {4!s} {5!s} --output {6!s}'.format(
           tool_path_check(TOOLS_DICT['transrate'].full_exe[0]), assembly_path, lefts,
           rights, cpu_cap, reference, transrate_dir)
    name = 'transrate_' + assembly_name
    out, err = GEN_LOGS(name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, cpu=cpu_cap, stdout=out, stderr=err)


def assembly_stats_task(out_dir,assembly,tasks):
    ''' Defines assembly_stats task. Uses PATH_DIR, PATH_SCRIPTS, NAME_ASSEMBLY.
        Params :
            tasks - a list of tasks that this task is dependant on (trinity_task)
    '''
    trgs = ['{0!s}/assembly_stats.json'.format(out_dir)]
    cmd = 'python {0!s}/assembly_stats.py {1!s} > {2!s}'.format(PATH_SCRIPTS,assembly,trgs[0])
    name = 'assembly_stats'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


##################___Annotation_Tasks___##################

def gene_trans_map_task(path_assembly,assembly_name,out_dir,tasks):
    '''    Defines gene_trans_map task. Uses NAME_ASSEMBLY, PATH_DIR, PATH_GENE_TRANS_MAP.
        Params :
            tasks - a list of tasks that this task is dependant on (trinity_task) 
    '''
    trgs = ['{0!s}/{1!s}.gene_trans_map'.format(out_dir, assembly_name)]
    cmd = '{0!s} {1!s} > {2!s}'.format(tool_path_check(TOOLS_DICT['trinity'].full_exe[1]),path_assembly,trgs[0])
    name = 'gene_trans_map'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def blastx_task(path_db, cpu_cap, tasks):
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
            'send evalue bitscore stitle slen" -evalue 0.0001 > {5!s}'
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
    cmd = ("cd {0!s}; {1!s} --transcriptome {2!s}/{3!s}.fasta  --path_to_rnammer {4!s} "
            "--org_type euk; cd -").format(GEN_PATH_ANNOTATION_FILES(),PATH_RNAMMER_PL,
            GEN_PATH_DIR(),NAME_ASSEMBLY,PATH_RNAMMER)
    name = 'rnammer'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def predict_orfs_task(cpu_cap, tasks):
    '''    defines predict_orfs task. Uses NAME_ASSEMBLY, PATH_DIR, transdecoder_tool.full_exe.
        Params : 
            cpu_cap - the number of threads to be used by task
            tasks - a list of tasks that this task is dependant on (trinity_task)
        blasp, pfamm, tmhmm, signalp are children

        *transdecoder.* are targets
    '''
    path_transdecoder_output = GEN_PATH_ANNOTATION_FILES()+'/transdecoder'
    trgs = ['{0!s}/{1!s}.fasta.transdecoder.pep'.format(path_transdecoder_output,NAME_ASSEMBLY)]
    cmd = ("mkdir -p {0!s}; cd {0!s}; {1!s} -t {2!s} --workdir {0!s} --CPU {3!s} 2>&1 | "
            "tee {0!s}/transDecoder_no_pfam_log").format(path_transdecoder_output,
            #tool_path_check(TOOLS_DICT['transdecoder'].full_exe[0]),GEN_PATH_ASSEMBLY(),cpu_cap)
            tool_path_check(PATH_TRANSDECODER),GEN_PATH_ASSEMBLY(),cpu_cap)
    name = 'TransDecoder'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def signalp_task(tasks):
    trgs = ['{0!s}/{1!s}.signalp'.format(GEN_PATH_ANNOTATION_FILES(),NAME_ASSEMBLY)]
    cmd = '{2!s} -f short -n {1!s} {0!s}'.format(GEN_PATH_PEP(),trgs[-1],PATH_SIGNALP)
    name = 'signalp'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def blastp_task(path_db, cpu_cap, tasks):
    '''    Defines a task for running blastp. Uses PATH_DIR, PATH_BLASTP.
        Params : 
            path_db - path to blastp databse
            cpu_cap - number of threads used by task
            tasks - a list of tasks that this task is dependant on (predict_orfs_task)
    '''
    db_name = os.path.basename(path_db).split('.')[0]
    trgs = ["{0!s}/{1!s}_{2!s}.blastp".format(GEN_PATH_ANNOTATION_FILES(), NAME_ASSEMBLY,db_name)]
    cmd = ('{0!s} -query {1!s} -db {2!s} -num_threads {3!s} -max_target_seqs 1 '
            '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart '
            'send evalue bitscore stitle slen" -evalue 0.001 > {4!s}'
            ).format(PATH_BLASTP,GEN_PATH_PEP(),path_db,cpu_cap,trgs[0])
    name = 'blastp_{0!s}'.format(db_name)
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def tmhmm_task(tasks):
    trgs = ['{0!s}/{1!s}.tmhmm'.format(GEN_PATH_ANNOTATION_FILES(),NAME_ASSEMBLY)]
    cmd = 'cd {0!s}; {1!s} --short < {2!s} > {3!s}'.format(GEN_PATH_ANNOTATION_FILES(),PATH_TMHMM,GEN_PATH_PEP(),trgs[0])
    name = 'tmhmm'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def pfam_task(cpu_cap, tasks):
    trgs = ['{0!s}/{1!s}.pfam'.format(GEN_PATH_ANNOTATION_FILES(),NAME_ASSEMBLY)]
    cmd = '{4!s} --cpu {0!s} --domtblout {1!s} {2!s} {3!s}'.format(cpu_cap,
            trgs[0],PATH_PFAM_DATABASE,GEN_PATH_PEP(),tool_path_check(TOOLS_DICT['hmmer'].full_exe[0]))
    name = 'pfam'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def annot_table_task(opts, tasks):
    suffixes = ['annotation.txt','annotation_by_gene.txt']
    trgs = ['{0!s}/{1!s}_{2!s}'.format(GEN_PATH_DIR(),NAME_ASSEMBLY,sufx) for sufx in suffixes]
    cmd = (
        'python {0!s}/annot_table_main.py --fasta {1!s} --outfile {2!s}/{3!s} '
        '--ko2path {4!s}/orthology_pathway.list --sp2enzyme '
        '{4!s}/swiss_enzyme.list --enzyme2path {4!s}/enzyme_pathway.list '
        '--pfam2enzyme {4!s}/pfam_enzyme.list --go2path {4!s}/go_pathway.txt '
        '--nog2function {4!s}/allKOG_functional_info.txt '
        '--go2slim {4!s}/goslim_generic.obo --sp2ko {4!s}/idmapping.KO '
        '--sp2nog {4!s}/idmapping.eggNOG --sp2ortho {4!s}/idmapping.orthodb '
        '--sp2bioc {4!s}/idmapping.biocyc --sp2goentrez '
        '{4!s}/idmapping_selected.tab ').format(
        PATH_SCRIPTS, GEN_PATH_ASSEMBLY(), GEN_PATH_DIR(), NAME_ASSEMBLY,
        PATH_DATABASES)
    cmd += ' '.join(['--'+k+' '+opts[k] for k in opts])
    name = 'build_annotation_table'
    out, err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def kegg_task(tasks):
    '''    Defines the kegg task. Uses PATH_SCRIPTS, 
        Params : 
    '''
    trgs = ['{0!s}/ko01100.pdf'.format(GEN_PATH_ANNOTATION_FILES()),
            '{0!s}/ko01100_KO.txt'.format(GEN_PATH_ANNOTATION_FILES())]
    cmd = ('python {0!s}/color_pathways2.py --path ko01100 --transcriptomeKO '
            '{1!s}/uniq_ko_annots.txt --output {2!s}').format(PATH_SCRIPTS,PATH_DATABASES,GEN_PATH_ANNOTATION_FILES())
    name='draw_kegg_maps'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def pipeplot_task(annotation_table, tasks):
    trgs = []
    cmd = 'mkdir -p {0!s}/plots ; cd {0!s}/plots ; python {1!s}/pipePlot.py -i {2!s} ;'.format(
            GEN_PATH_ANNOTATION_FILES(),PATH_SCRIPTS,annotation_table)
    name = 'pipeplot'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def diamondX_task(ref, cpu_cap, tasks):
    base_ref = os.path.basename(ref)
    #trgs = ['{0!s}/diamond_{1!s}.blastx'.format(GEN_PATH_ANNOTATION_FILES(), base_ref)]
    trgs = ['{0!s}/{1!s}_x_{2!s}.diamond_blastx'.format(GEN_PATH_ANNOTATION_FILES(), NAME_ASSEMBLY, base_ref)]
    pseudo_trgs = ['{0!s}/diamond_{1!s}_blastx'.format(GEN_PATH_ANNOTATION_FILES(), base_ref)]
    cmd = ('{0!s} blastx --db {1!s} --query {2!s} --daa {3!s} --tmpdir {4!s} '
           '--max-target-seqs 20 --sensitive --threads {5!s} --evalue 0.001; {0!s} view '
           '--daa {3!s}.daa --out {6!s};').format(
           tool_path_check(TOOLS_DICT['diamond'].full_exe[0]), ref, GEN_PATH_ASSEMBLY(), pseudo_trgs[0], GEN_PATH_ANNOTATION_FILES(),
           cpu_cap, trgs[0])
    name = 'diamond_blastx_'+base_ref
    out, err = GEN_LOGS(name)
    return Task(command=cmd, dependencies=tasks, cpu=cpu_cap, targets=trgs, name=name, stdout=out, stderr=err)


def diamondP_task(ref, cpu_cap, tasks):
    base_ref = os.path.basename(ref)
    #trgs = ['{0!s}/diamond_{1!s}.blastp'.format(GEN_PATH_ANNOTATION_FILES(), base_ref)]
    trgs = ['{0!s}/{1!s}_x_{2!s}.diamond_blastp'.format(GEN_PATH_ANNOTATION_FILES(), NAME_ASSEMBLY, base_ref)]
    pseudo_trgs = ['{0!s}/diamond_{1!s}_blastp'.format(GEN_PATH_ANNOTATION_FILES(), base_ref)]
    cmd = ('{0!s} blastp --db {1!s} --query {2!s} --daa {3!s} --tmpdir {4!s} '
           '--max-target-seqs 20 --sensitive --threads {5!s} --evalue 0.001; {0!s} view '
           '--daa {3!s}.daa --out {6!s};').format(
           tool_path_check(TOOLS_DICT['diamond'].full_exe[0]), ref, GEN_PATH_PEP(), pseudo_trgs[0], GEN_PATH_ANNOTATION_FILES(),
           cpu_cap, trgs[0])
    name = 'diamond_blastp_'+base_ref
    out, err = GEN_LOGS(name)
    return Task(command=cmd, dependencies=tasks, cpu=cpu_cap, targets=trgs, name=name, stdout=out, stderr=err)


def blast_augment_task(db, blast, tasks):
    id2name = db+'.stitle'
    trgs = ['{0!s}_ex'.format(blast)]
    cmd = 'python {0!s}/addStitleToBlastTab.py --db2Name {1!s} --blast {2!s} > {3!s}'.format(
           PATH_SCRIPTS, id2name, blast, trgs[0])
    name = 'Blast_Augmentation_'+os.path.basename(blast)
    out, err = GEN_LOGS(name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def build_bowtie_task(path_assembly, assembly_name, out_dir, tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.1.bt2'.format(out_dir,assembly_name)]  
    cmd = 'bowtie2-build --offrate 1 -f {1!s} {2!s}/{3!s}'.format(
            PATH_BOWTIE2, path_assembly, out_dir, assembly_name) 
    name = 'build_bowtie'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def bowtie2_unpaired_task(bowtie2_index,out_dir,fastq,out_name,opt,cpu_cap,tasks):
    '''
    '''
    opts = ['-a -t --end-to-end', '-t --local']
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir,out_name)]
    cmd = ('bowtie2 {1!s} -L {2!s} -N 1 --threads {3!s} -x {4!s} -U '
            '{5!s} | samtools view -Sb - > {6!s} ').format(PATH_BOWTIE2,
            opts[opt],22,cpu_cap,bowtie2_index,fastq,trgs[0])
    name = 'bowtie2_'+out_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def bowtie2_task(bowtie2_index,out_dir,fastq1,fastq2,out_name,opt,cpu_cap,tasks):
    '''    
    '''
    opts = ['-a -t --end-to-end', '-t --local']
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir,out_name)]
    cmd = ('bowtie2 {1!s} -L {2!s} -N 1 --maxins 800 --threads {3!s} -x {4!s} -1 '
            '{5!s} -2 {6!s} | samtools view -Sb - > {7!s} ').format(PATH_BOWTIE2,
            opts[opt],22,cpu_cap,bowtie2_index,fastq1,fastq2,trgs[0])
    name = 'bowtie2_'+out_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def express_task(assembly_path,out_dir,out_name,bam_input,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.xprs'.format(out_dir,out_name)]
    cmd = ('mkdir {1!s}/{2!s}; {0!s} --output-dir {1!s}/{2!s} {3!s} {4!s}; mv '
            '{1!s}/{2!s}/results.xprs {5!s}; rm -rf {1!s}/{2!s};').format(
            PATH_EXPRESS,out_dir,out_name,assembly_path,bam_input,trgs[0])
    name = 'express_'+out_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def counts_to_table_task(gene_trans_map,out_dir,express_files,out_name,flag,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.countsTable'.format(out_dir,out_name)]
    count_str = ' '.join(['--counts {0!s}'.format(f) for f in express_files])
    cmd = ('python {0!s}/counts_to_table2.py --out {1!s} --inDir {2!s} '
            '--outDir {2!s} {3!s} {4!s} --geneTransMap {5!s}').format( 
            PATH_SCRIPTS,out_name,out_dir,flag,count_str,gene_trans_map)
    name = 'counts_to_table_'+flag
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def sam_sort_task(out_dir,bam_file,out_name,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir,out_name)]
    cmd = 'samtools sort {0!s} {1!s}/{2!s}'.format(bam_file,out_dir,out_name)
    name = 'sam_sort_'+out_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def assembly_to_bed_task(path_assembly,assembly_name,out_dir, tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.bed'.format(out_dir,assembly_name)]
    cmd = 'python {0!s}/fasta_to_bed_count_length.py {1!s} {2!s}'.format(
            PATH_SCRIPTS,path_assembly,trgs[0])
    name = 'fasta_to_bed'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def intersect_bed_task(out_dir,bam_file,bed_reference,output_name,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.bed'.format(out_dir,output_name)]
    cmd = '{0!s} intersect -abam {1!s} -b {2!s} -wb -bed > {3!s}'.format(PATH_BEDTOOLS,bam_file,bed_reference,trgs[0])
    name = 'intersect_bed_'+output_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def deseq2_task(out_dir,counts_to_table_results,sample_info,basename,model,tasks):
    '''
    '''
    pseudo_model_temp = ''.join([c if c!=' ' else '_' for c in model])
    trgs = ['{0!s}/deseq2_{1!s}_{2!s}/'.format(out_dir,basename,pseudo_model_temp)]
    cmd = 'Rscript {5!s}/deseq2.r --args {0!s} {1!s} {2!s} {3!s} {4!s}'.format(
            counts_to_table_results,sample_info,out_dir,basename,model,PATH_SCRIPTS)
    name = 'de_'+basename
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def build_salmon_task(path_assembly,assembly_name,out_dir,cpu_cap,tasks):
    trgs = ['{0!s}/{1!s}_salmon'.format(out_dir, assembly_name)] 
    cmd = '{0!s} index --transcripts {1!s} --index {2!s}/{3!s}_salmon --threads {4!s} --type quasi'.format(tool_path_check(TOOLS_DICT['salmon'].full_exe[0]),path_assembly, out_dir, assembly_name, cpu_cap)
    name = 'build_salmon_' + assembly_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err, cpu=cpu_cap)


def salmon_gene_map_task(out_dir,assembly_name,gene_trans_map,tasks):
    ''' salmon requires gene_trans_map in reverse column order (transcript \t gene \n)'''
    trgs = ['{0!s}/{1!s}.trans_gene_map'.format(out_dir,assembly_name)] 
    cmd = 'join -t, -o 1.2,1.1 {0!s} {0!s} > {1!s}'.format(gene_trans_map, trgs[0]) 
    name = 'trans_gene_map'
    name = 'salmon_gene_map_task_' + assembly_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def salmon_task(index,left,right,out_name,gene_map,out_dir,cpu_cap,tasks):
    trgs = ['{0!s}/{1!s}/quant.sf'.format(out_dir,out_name)]
    cmd = '{0!s} quant -i {1!s} -l IU -1 {2!s} -2 {3!s} -o {4!s}/{5!s} --geneMap {6!s} -p {7!s} --extraSensitive; cp ' \
        '{4!s}/{5!s}/quant.sf {4!s}/{5!s}_quant.sf; cp {4!s}/{5!s}/quant.genes.sf {4!s}/{5!s}_quant.genes.sf'.format(
	tool_path_check(TOOLS_DICT['salmon'].full_exe[0]),index,left,right,out_dir,out_name,gene_map,cpu_cap)
    name = 'salmon_' + os.path.basename(index)
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def salmon_unpaired_task(index,unpaired,out_name,gene_map,out_dir,cpu_cap,tasks):
    trgs = []
    cmd = '{0!s} quant -i {1!s} -l U -r {2!s} -o {3!s}/{4!s} --geneMap {5!s} -p {6!s} --extraSensitive'.format(
            tool_path_check(TOOLS_DICT['salmon'].full_exe[0]),index,unpaired,out_dir,out_name,gene_map,cpu_cap)
    name = 'salmon_unpaired_' + os.path.basename(index)
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def build_kallisto_task(assembly_path, assembly_name,out_dir,tasks):
    trgs = []
    cmd = '{0!s} index -i {1!s}/{2!s}_kallisto {3!s}'.format(
            PATH_KALLISTO,out_dir,assembly_name,assembly_path)
    name = 'build_kallisto'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def kallisto_task(index,out_dir,out_name,left,right,tasks):
    trgs = []
    cmd = '{0!s} quant -i {1!s} -o {2!s}/{3!s} {4!s} {5!s}'.format(
            PATH_KALLISTO,index,out_dir,out_name,left,right)
    name = 'kallisto'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def build_blast_task(fasta,out_path,dbtype,tasks,log_flag=True):
    trgs = []
    #title doesn't seem to change the out name .. it's still xx.gz.psq, etc? CHECK.
    title = os.path.basename(out_path).split('.')[0]
    cmd = 'gunzip -c {0!s} | makeblastdb -in - -dbtype {2!s} -title {3!s} -out {1!s}'.format(fasta,out_path,dbtype,title)
    name = 'build_blast_'+title
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def build_diamond_task(fasta,out_path,tasks,log_flag=True):
    title = os.path.basename(out_path)
    trgs = ['{0!s}'.format(out_path + '.dmnd')] 
    cmd = '{0!s} makedb --in {1!s} --db {2!s}'.format(TOOLS_DICT['diamond'].full_exe[0], fasta, out_path)
    name = 'build_diamond_'+ title
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)
 

def split_mito_task(blast_mt,tasks):
    trgs = ['{0!s}/mtDNA_contigs.fasta'.format(GEN_PATH_ASSEMBLY_FILES()),'{0!s}/no_mtDNA_contigs.fasta'.format(GEN_PATH_ASSEMBLY_FILES())]
    cmd = '{0!s}/split_fasta.py {1!s} {3!s} {2!s}/mtDNA_contigs.fasta {2!s}/no_mtDNA_contigs.fasta'.format(PATH_SCRIPTS,GEN_PATH_ASSEMBLY(),GEN_PATH_ASSEMBLY_FILES(),blast_mt)
    name = 'split_mito'
    out, err = GEN_LOGS(name) 
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def pfam_build_task(source, tasks, log_flag=True):
    trgs = [PATH_PFAM_DATABASE+'.h3f']
    cmd = 'cd {0!s} ; {1!s} -f {2!s};'.format(PATH_DATABASES, tool_path_check(TOOLS_DICT['hmmer'].full_exe[1]), source)
    name = 'hmmpress'
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def db2stitle_task(db, tasks, log_flag=True):
    base_db = os.path.basename(db)
    trgs = ['{0!s}/{1!s}.stitle'.format(PATH_DATABASES, base_db)]
    cmd = 'python {0!s}/fastaID2names.py --fasta {1!s} > {2!s}'.format(
           PATH_SCRIPTS, db, trgs[0])
    name = 'db2stitle_'+base_db
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)

def manage_db_task(fresh, nr_flag, uniref90_flag, busco_flags, blastplus_flag, cpu_cap, tasks, log_flag=True):
    trgs = [PATH_DATABASES]
    cmd = 'python {0!s}/manage_database.py'.format(PATH_SCRIPTS)
    if(fresh):
        cmd +=' --hard'
    if(nr_flag):
        cmd+=' --nr'
    if(uniref90_flag):
        cmd+=' --uniref90'
    if(len(busco_flags) != 0):
        cmd+=' --buscos '
        cmd+=','.join(busco_flags)
    if blastplus_flag:
        cmd+= ' --buildBlastPlus'
    name = 'db_manage'
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap)

def manage_tools_task(install, fresh, cpu_cap, tool_list, tasks, log_flag=True):
    trgs = [PATH_TOOLS]
    cmd = 'python {0!s}/manage_tools.py'.format(PATH_SCRIPTS)
    if(install):
        cmd+=' --install'
    if(fresh):
        cmd +=' --hard'
    cmd+= ' --tool ' + ' --tool '.join(tool_list)
    name = 'tools_manage'
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap)

def filter_task(assembly_path, assembly_name, out_dir, quant_file_list, tpm_threshold, tpm_column_index, tasks, log_flag=True):
    # TPM column index: transrate uses older salmon; use index =2. Newer salmon: index=3
    trgs = ['{0!s}/{1!s}_{2!s}tpm.fasta'.format(out_dir,assembly_name,tpm_threshold)]
    quants = ''.join(' --quant_files '+ x for x in quant_file_list) 
    cmd = 'python {0!s}/filter_contigs_by_tpm.py --assembly {1!s} --tpm {2!s} {3!s} --out {4!s} --tpm_column_index {5!s}'.format(PATH_SCRIPTS, assembly_path,tpm_threshold, quants, trgs[0], tpm_column_index)
    name = 'filt_{0!s}_{1!s}tpm'.format(assembly_name, tpm_threshold)
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)

