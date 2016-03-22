'''
'''

from tasks_v2 import Task
import os
from os.path import join, exists
import sys
import functions_general as fg
from external_tools import TOOLS_DICT
import re
import warnings

''' static db variables '''
PATH_BUSCO_REFERENCE = join(fg.PATH_DATABASES, 'busco')

def cegma_task(out_dir,assembly,cpu_cap, tasks):
    '''    Defines the cegma task. Uses PATH_DIR, PATH_CEGMA, NAME_ASSEMBLY.
        Params :
            cpu_cap - number of threads to be used by cegma
            tasks - a list of tasks that this task is dependant on (trinity_task)
    '''
    assembly_name = os.path.basename(assembly).split('.fa')[0]
    trgs = ['{0!s}/{1!s}.completeness_report'.format(out_dir,assembly_name)]
    cmd = '{0!s} -g {1!s} -v -o {3!s}/{2!s} -T {4!s}'.format(fg.tool_path_check(TOOLS_DICT['cegma'].full_exe[0]),
            assembly,assembly_name,out_dir,cpu_cap)
    name = 'cegma'
    out,err = fg.GEN_LOGS(name)
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
            ).format(out_dir,fg.tool_path_check(TOOLS_DICT['busco_plant'].full_exe[0]),reference_name,assembly_name,assembly_path,
            PATH_BUSCO_REFERENCE,cpu_cap)
    name = 'busco_'+ reference_name + '_' + assembly_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,cpu=cpu_cap,stdout=out,stderr=err)


def transrate_dep_generator(reads_dir,transrate_task, lefts, rights, reference, assembly_path, cpu_cap, transrate_dir, other_dependencies):

    def ret():
        for t in other_dependencies:
            try:
                if( not d.finished()):
                    return False
            except Task.ExitCodeException:
                return False
        assembly_files = sorted(os.listdir(reads_dir))
        assembly_files = [os.path.join(reads_dir, f) for f in assembly_files]
        new_lefts = [[g for g in assembly_files if(os.path.basename(f) in g)] for f in lefts]
        new_lefts = [k[0] for k in new_lefts if(len(k) > 0)]
        new_rights = [[g for g in assembly_files if(os.path.basename(f) in g)] for f in rights]
        new_rights = [k[0] for k in new_rights if(len(k) > 0)]
#        new_singles = [[g for g in assembly_files if(os.path.basename(f) in g)] for f in singles]
#        new_singles = [k[0] for k in new_singles if(len(k) > 0)]
        if(len(new_lefts) == len(lefts) and len(new_rights) == len(rights)): #and len(new_singles) == len(singles)
#            and len(new_lefts)+len(new_singles) != 0):
            new_lefts = ','.join(new_lefts)#+new_singles)
            new_rights = ','.join(new_rights) 
            new_lefts = '--left '+new_lefts if(len(new_lefts) > 0) else ''
            new_rights = '--right '+new_rights if(len(new_rights) > 0) else ''
            cmd = '{0!s} --assembly {1!s} {2!s} {3!s} --threads {4!s} {5!s} --output {6!s}'.format(
                   fg.tool_path_check(TOOLS_DICT['transrate'].full_exe[0]), assembly_path, new_lefts,
                   new_rights, cpu_cap, reference, transrate_dir)
            transrate_task.command = cmd
        else:
            warnings.warn('Unable to match input files with trimmed output. Continuing transrate using input files instead.')
        return True
    return ret


def transrate_task(reads_dir, assembly_path, assembly_name,lefts, rights, out_dir, transrate_dir, cpu_cap, tasks, reference = ''): #, cpu_cap, tasks):
    trgs = ['{0!s}/assemblies.csv'.format(transrate_dir),'{0!s}/{1!s}/good.{1!s}.fasta'.format(transrate_dir,assembly_name),'{0!s}/{1!s}/{1!s}.fasta_quant.sf'.format(transrate_dir,assembly_name)]
    orig_lefts = lefts
    orig_rights = rights
    lefts = ','.join(lefts)
    rights = ','.join(rights) 
    lefts = '--left '+lefts if(len(lefts) > 0) else ''
    rights = '--right '+rights if(len(rights) > 0) else ''
    reference = '--reference ' + reference if(reference != '') else ''
    cmd = '{0!s} --assembly {1!s} {2!s} {3!s} --threads {4!s} {5!s} --output {6!s}'.format(
           fg.tool_path_check(TOOLS_DICT['transrate'].full_exe[0]), assembly_path, lefts,
           rights, cpu_cap, reference, transrate_dir)
    name = 'transrate_' + assembly_name
    out, err = fg.GEN_LOGS(name)
    temp_task = Task(command=cmd, dependencies=[], targets=trgs, name=name, cpu=cpu_cap, stdout=out, stderr=err, max_wall_time=720)
    deps = transrate_dep_generator(reads_dir, temp_task, orig_lefts, orig_rights, reference, assembly_path, cpu_cap, transrate_dir, tasks)
    temp_task.dependencies = [deps]
    return temp_task


def assembly_stats_task(out_dir,assembly,tasks):
    ''' Defines assembly_stats task. Uses PATH_DIR, PATH_SCRIPTS, NAME_ASSEMBLY.
        Params :
            tasks - a list of tasks that this task is dependant on (trinity_task)
    '''
    trgs = ['{0!s}/assembly_stats.json'.format(out_dir)]
    cmd = 'python {0!s}/assembly_stats.py {1!s} > {2!s}'.format(fg.PATH_SCRIPTS,assembly,trgs[0])
    name = 'assembly_stats_' + os.path.basename(assembly)
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def filter_task(assembly_path, assembly_name, out_dir, quant_file_list, tpm_threshold, tpm_column_index, tasks, log_flag=True):
    # TPM column index: transrate uses older salmon; use index =2. Newer salmon: index=3
    trgs = ['{0!s}/{1!s}_{2!s}tpm.fasta'.format(out_dir,assembly_name,tpm_threshold)]
    quants = ''.join(' --quant_files '+ x for x in quant_file_list) 
    cmd = 'python {0!s}/filter_contigs_by_tpm.py --assembly {1!s} --tpm {2!s} {3!s} --out {4!s} --tpm_column_index {5!s}'.format(fg.PATH_SCRIPTS, assembly_path,tpm_threshold, quants, trgs[0], tpm_column_index)
    name = 'filt_{0!s}_{1!s}tpm'.format(assembly_name, tpm_threshold)
    out, err = fg.GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)

