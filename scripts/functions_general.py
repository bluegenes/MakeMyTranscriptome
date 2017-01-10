import mmt_defaults as statics
from tasks_v2 import Task
import os
import sys
import warnings
if(sys.version[0] == '3'):
    from shutil import which
else:
    from py2_which import which_python2 as which


"""
def stufF:
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


    class Output_Path_Vars:

        def __init__(self, basename, out_dir=statics.PATH_ASSEMBLIES):
            join = os.path.join
            self.path_dir = join(out_dir, basename)
            self.path_assembly_files = join(self.path_dir, 'assembly_files')
            self.path_quality_files = join(self.path_dir, 'quality_files')
            self.path_annotation_files = join(self.path_dir, 'annotation_files')
            self.path_expression_files = join(self.path_dir, 'expression_files')
            self.path_filter_files = join(self.path_dir, 'filtered_assemblies')
            self.path_logs = join(self.path_dir, 'log_files')
            self.path_assembly = join(self.path_dir, basename+'.fasta')
            self.path_gene_trans_map = join(self.path_assembly_files, basename+'.gene_trans_map')
            self.path_transdecoder_dir = join(self.path_annotation_files, 'transdecoder')
            self.path_transrate_dir = join(self.path_quality_files, 'transrate')
            self.path_pep = join(self.path_transdecoder_dir, basename+'.fasta.transdecoder.txt')
            self.path_annot_table = join(self.path_dir, basename+'_annotation.txt')

"""

# opc is the output_path_class object


def gen_logs(log_dir, name):
    stderr = os.path.join(log_dir, name + '.err_log')
    stdout = os.path.join(log_dir, name + '.out_log')
    return (stdout, stderr)


def round_div(cpu, k):
    ret = int(round(float(cpu) / k))
    return ret


def tool_path_check(full_exe):
    name = os.path.basename(full_exe)
    if(os.path.exists(full_exe)):
        return full_exe
    elif(which(name)):
        err_string = ('MMT has not installed {0!s}. MMT will use the version '
                      'found in your path variable instead.').format(name)
        warnings.warn(err_string)
        return name
    else:
        err_string = ('INSTALLATION ERROR : MMT has not installed {0!s} and '
                      'did not find this program in your path variable. Steps '
                      'requiring {0!s} will fail.').format(name)
        warnings.warn(err_string)
        return name


def make_dir_task(path, tasks=[]):
    trgs = [path]
    cmd = 'mkdir -p {0!s}'.format(path)
    name = 'mdkir_' + os.path.basename(path)
    return Task(command=cmd, dependencies=[t for t in tasks], targets=trgs, stdout=os.devnull, stderr=os.devnull, name=name)


def cp_task(source, target, tasks):
    trgs = [target]
    cmd = "cp {0!s} {1!s}".format(source, target)
    name = 'cp_{0!s}_{1!s}'.format(os.path.basename(source), os.path.basename(target))
    return Task(command=cmd, dependencies=tasks, targets=trgs, stdout=os.devnull, stderr=os.devnull, name=name)


def build_dir_task(opc, tasks):
    trgs = [opc.path_dir, opc.path_assembly_files, opc.path_quality_files,
            opc.path_annotation_files, opc.path_filter_files,
            opc.path_expression_files, opc.path_logs]
    cmd = 'mkdir -p {0!s}'.format(' '.join(trgs))
    return Task(command=cmd, dependencies=tasks, targets=trgs, stdout=os.devnull, stderr=os.devnull)


def cp_assembly_task(path_assembly, source, tasks):
    # may fail if the cp fails due to the target already existing
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    trgs = ['{0!s}'.format(path_assembly)]
    cmd = 'cp {0!s} {1!s}'.format(source, trgs[0])
    name = 'setting_fasta_' + assembly_name
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name)


def split_mito_task(opc, blast_mt, tasks):
    trgs = ['{0!s}/mtDNA_contigs.fasta'.format(opc.path_assembly_files),
            '{0!s}/no_mtDNA_contigs.fasta'.format(opc.path_assembly_files)]
    cmd = ('{0!s}/split_fasta.py {1!s} {3!s} {2!s}/mtDNA_contigs.fasta '
           '{2!s}/no_mtDNA_contigs.fasta').format(
           statics.PATH_UTIL, opc.path_assembly, opc.path_assembly_files, blast_mt)
    name = 'split_mito'
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


'''
def db2stitle_task(db, tasks, log_flag=True, opc=None):
    # does this belong here This seems like a database function??
    base_db = os.path.basename(db)
    trgs = ['{0!s}/{1!s}.stitle'.format(statics.PATH_DATABASES, base_db)]
    cmd = 'python {0!s}/fastaID2names.py --fasta {1!s} > {2!s}'.format(
           statics.PATH_UTIL, db, trgs[0])
    name = 'db2stitle_'+base_db
    out, err = gen_logs(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)
'''


def manage_tools_task(install, fresh, cpu_cap, tool_list, tasks, log_flag=True, opc=None):
    trgs = [statics.PATH_TOOLS]
    cmd = 'python {0!s}/manage_tools.py'.format(statics.PATH_SCRIPTS)
    if(install):
        cmd += ' --install'
    if(fresh):
        cmd += ' --hard'
    cmd += ' --tool ' + ' --tool '.join(tool_list)
    name = 'tools_manage'
    base_log = os.path.join(statics.PATH_TOOLS, name)
    out, err = base_log + '.out', base_log + '.err'
    out, err = (out, err) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap)

def filter_task(assembly_path, assembly_name, out_dir, quant_file_list, tpm_threshold, tpm_column_index, tasks, log_flag=True, opc=None):
    '''TPM column index: transrate uses older salmon; use index =2. Newer salmon: index=3'''
    trgs = ['{0!s}/{1!s}_{2!s}tpm.fasta'.format(out_dir, assembly_name, tpm_threshold)]
    quants = ''.join(' --quant_files ' + x for x in quant_file_list)
    cmd = ('python {0!s}/filter_contigs_by_tpm.py --assembly {1!s} --tpm {2!s} {3!s} '
           '--out {4!s} --tpm_column_index {5!s}').format(
           statics.PATH_UTIL, assembly_path, tpm_threshold, quants, trgs[0],
           tpm_column_index)
    name = 'filt_{0!s}_{1!s}tpm'.format(assembly_name, tpm_threshold)
    out, err = os.devnull, os.devnull
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)
