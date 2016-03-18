'''
'''
from tasks_v2 import Task
import os
from os.path import join, exists, dirname, abspath
import sys
if(sys.version[0] == '3'):
    from shutil import which
else:
    from py2_which import which_python2 as which

import re
import warnings

''' name variables '''
NAME_ASSEMBLY = 'myassembly'
NAME_OUT_DIR = 'mmt_test_output'

''' static path variables '''
PATH_ROOT = dirname(dirname(abspath(__file__)))
PATH_TEST = join(PATH_ROOT, 'test_data')
PATH_SCRIPTS = join(PATH_ROOT, 'scripts')
PATH_DATABASES = join(PATH_ROOT, 'databases')
PATH_ASSEMBLIES = join(PATH_ROOT, 'assemblies')
PATH_TOOLS = join(PATH_ROOT, 'external_tools')

''' static db variables '''
PATH_BUSCO_REFERENCE = join(PATH_DATABASES, 'busco')
PATH_PFAM_DATABASE = '{0!s}/pfam/Pfam-A.hmm'.format(PATH_DATABASES)
PATH_NR = join(PATH_DATABASES, 'nr', 'nr')
PATH_SWISS_PROT = join(PATH_DATABASES, 'uniprot_sprot', 'uniprot_sprot')
PATH_UNIREF90 = join(PATH_DATABASES, 'uniref90', 'uniref90')
PATH_NOG_CATEGORIES = join(PATH_DATABASES, 'nog_categories')


# Dynamic path variable functions
def GEN_PATH_DIR(): return os.path.join(PATH_ASSEMBLIES, NAME_OUT_DIR)

def GEN_PATH_ASSEMBLY_FILES(): return os.path.join(GEN_PATH_DIR(), 'assembly_files')

def GEN_PATH_QUALITY_FILES(): return os.path.join(GEN_PATH_DIR(), 'quality_files')

def GEN_PATH_ANNOTATION_FILES(): return os.path.join(GEN_PATH_DIR(), 'annotation_files')

def GEN_PATH_EXPRESSION_FILES(): return os.path.join(GEN_PATH_DIR(), 'expression_files')

def GEN_PATH_FILTER_FILES(): return os.path.join(GEN_PATH_DIR(), 'filtered_assemblies')

def GEN_PATH_LOGS(): return os.path.join(GEN_PATH_DIR(), 'log_files')

def GEN_PATH_ASSEMBLY(): return os.path.join(GEN_PATH_DIR(), NAME_ASSEMBLY+'.fasta')

def GEN_PATH_GENE_TRANS_MAP(): return os.path.join(GEN_PATH_ASSEMBLY_FILES(), NAME_ASSEMBLY+'.gene_trans_map')

def GEN_PATH_TRANSDECODER_DIR(): return os.path.join(GEN_PATH_ANNOTATION_FILES(), 'transdecoder')

def GEN_PATH_TRANSRATE_DIR(): return os.path.join(GEN_PATH_QUALITY_FILES(), 'transrate')

def GEN_PATH_PEP(): return os.path.join(GEN_PATH_TRANSDECODER_DIR(), NAME_ASSEMBLY+'.fasta.transdecoder.pep')

def GEN_PATH_ANNOT_TABLE(): return os.path.join(GEN_PATH_DIR(), NAME_ASSEMBLY+'annotation.txt')



def GEN_LOGS(x): return (os.path.join(GEN_PATH_LOGS(), x+'.out_log'),
                         os.path.join(GEN_PATH_LOGS(), x+'.err_log'))

def tool_path_check(full_exe):
    name = os.path.basename(full_exe)
    if(os.path.exists(full_exe)):
        return full_exe
    elif(which(name)):
        warnings.warn('MMT has not installed '+name+'. MMT will use the version found in your path variable instead.')
        return name
    else:
        warnings.warn('INSTALLATION ERROR : MMT has not installed '+name+' and did not find this program in your path variable. Steps requiring '+name+' will not be performed.') 


def build_dir_task(tasks):
    '''
    '''
    trgs = [GEN_PATH_DIR(), GEN_PATH_ASSEMBLY_FILES(), GEN_PATH_QUALITY_FILES(), GEN_PATH_ANNOTATION_FILES(),
            GEN_PATH_FILTER_FILES(), GEN_PATH_EXPRESSION_FILES(), GEN_PATH_LOGS()]
    cmd = ' '.join(['mkdir -p {0!s};'.format(d) for d in trgs])
    return Task(command=cmd,dependencies=tasks,targets=trgs,stdout=os.devnull,stderr=os.devnull)


def cp_assembly_task(path_assembly, source, tasks):
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    '''    Defines task used to initialize an assembly when running on fasta
        files. Uses GEN_PATH_DIR() and NAME_ASSEMBLY.
        Params :
            source - The path to the source fasta that should be used for analsis
            tasks - a list of tasks that this task is dependant on.
    '''
    trgs = ['{0!s}'.format(path_assembly)]
    cmd = 'cp {0!s} {1!s}'.format(source, trgs[0]) 
    name = 'setting_fasta_' + assembly_name
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name)


def split_mito_task(blast_mt,tasks):
    trgs = ['{0!s}/mtDNA_contigs.fasta'.format(GEN_PATH_ASSEMBLY_FILES()),'{0!s}/no_mtDNA_contigs.fasta'.format(GEN_PATH_ASSEMBLY_FILES())]
    cmd = '{0!s}/split_fasta.py {1!s} {3!s} {2!s}/mtDNA_contigs.fasta {2!s}/no_mtDNA_contigs.fasta'.format(PATH_SCRIPTS,GEN_PATH_ASSEMBLY(),GEN_PATH_ASSEMBLY_FILES(),blast_mt)
    name = 'split_mito'
    out, err = GEN_LOGS(name) 
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


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
    '''TPM column index: transrate uses older salmon; use index =2. Newer salmon: index=3'''
    trgs = ['{0!s}/{1!s}_{2!s}tpm.fasta'.format(out_dir,assembly_name,tpm_threshold)]
    quants = ''.join(' --quant_files '+ x for x in quant_file_list) 
    cmd = 'python {0!s}/filter_contigs_by_tpm.py --assembly {1!s} --tpm {2!s} {3!s} --out {4!s} --tpm_column_index {5!s}'.format(PATH_SCRIPTS, assembly_path,tpm_threshold, quants, trgs[0], tpm_column_index)
    name = 'filt_{0!s}_{1!s}tpm'.format(assembly_name, tpm_threshold)
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)

