'''
'''

from tasks_v2 import Task
import os
from external_tools import TOOLS_DICT
from functions_general import tool_path_check
import mmt_defaults as statics


''' static db variables
PATH_PFAM_DATABASE = '{0!s}/pfam/Pfam-A.hmm'.format(fg.PATH_DATABASES)
PATH_NR = join(fg.PATH_DATABASES, 'nr', 'nr')
PATH_SWISS_PROT = join(fg.PATH_DATABASES, 'uniprot_sprot', 'uniprot_sprot')
PATH_UNIREF90 = join(fg.PATH_DATABASES, 'uniref90', 'uniref90')
PATH_NOG_CATEGORIES = join(fg.PATH_DATABASES, 'nog_categories')
'''


def gen_db_logs(name):
    base_log_name = os.path.join(statics.PATH_DATABASE_LOGS, name)
    out, err = base_log_name+'.out', base_log_name+'.err'
    return out, err


def db2stitle_task(db, tasks):
    base_db = os.path.basename(db)
    trgs = ['{0!s}/{1!s}.stitle'.format(statics.PATH_DATABASES, base_db)]
    cmd = 'python {0!s}/fastaID2names.py --fasta {1!s} > {2!s}'.format(
           statics.PATH_SCRIPTS, db, trgs[0])
    name = 'db2stitle_' + base_db
    out, err = gen_db_logs(name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def build_blast_task(path_db, out_dir, dbtype, tasks):
    trgs = []
    # title doesn't seem to change the out name .. it's still xx.gz.psq, etc? CHECK.
    title = os.path.basename(path_db).split('.')[0]
    cmd = 'gunzip -c {0!s} | {1!s} -in - -dbtype {2!s} -title {3!s} -out {4!s}'.format(
          path_db, tool_path_check(TOOLS_DICT['blast'].full_exe[0]), dbtype, title, out_dir)
    name = 'build_blastplus_db_' + title
    out, err = gen_db_logs(name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def build_diamond_task(path_db_fasta, out_path, tasks):
    ''' Is there a reason that we aren't checking for the installation of daimond?
    '''
    title = os.path.basename(out_path)
    trgs = ['{0!s}'.format(out_path + '.dmnd')]
    cmd = '{0!s} makedb --in {1!s} --db {2!s}'.format(
          tool_path_check(TOOLS_DICT['diamond'].full_exe[0]), path_db_fasta, out_path)
    name = 'build_diamond_' + title
    out, err = gen_db_logs(name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def pfam_build_task(source, out_root_path, tasks):
    ''' Trgs seem to be declared without respect to input.
    '''
    trgs = [os.path.join(statics.PATH_PFAM_DIR, os.path.basename(source))+'.h3f']
    cmd = 'cd {0!s} ; {1!s} -f {2!s};'.format(
          statics.PATH_DATABASES, tool_path_check(TOOLS_DICT['hmmer'].full_exe[1]),
          source)
    name = 'hmmpress_' + os.path.basename(source)
    out, err = gen_db_logs(name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def download_task(url, install_location, ftype, tasks):
    trgs = [install_location]
    cmd = 'python {0!s}/url_retrieve.py {1!s} --target {2!s} --type {3!s}'.format(
           statics.PATH_SCRIPTS, url, install_location, ftype)
    name = 'download_{0!s}'.format(os.path.basename(install_location))
    out, err = gen_db_logs(name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)
