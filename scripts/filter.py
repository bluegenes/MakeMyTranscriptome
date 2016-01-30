import argparse
import os
from os.path import join,exists
import sys
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf
import time


# dependency_set = trinity. but also transrate, for the transrate_good filtering...
def gen_filter_supervisor(fq1, fq2, unpaired,  dependency_set, busco_refs, cpu=12, cegma_flag=False, transrate_ref='', filter_tpm=0.5):
    tasks = []
    transrate_good_assembly = join(tf.GEN_PATH_TRANSRATE_DIR(),NAME_ASSEMBLY,'good.'+NAME_ASSEMBLY+'.fasta')
    tr_good_build_salmon = tf.build_salmon_task(cpu,[])
    salmon_transrate_good = tf.salmon_task()
    salmon_full = tf.salmon_task()

    filter_transrate_good = tf.filter_task(transrate_good_assembly, filter_tpm, []) # transrate is a dependency - can it depend on the file existing?
    tasks.append(filter_transrate_good)
    filter_full_assembly = tf.filter_task(tf.GEN_PATH_ASSEMBLY(), filter_tpm, [])  #trinity is dependency .. can it depend on the file?
    tasks.append(filter_full_assembly)
    assemblies = [transrate_good_assembly, filter_transrate_good_task.targets[0], filter_full_assembly.targets[0]]  
    for assemb in assemblies:
        tasks.append(tf.assembly_stats_task(assemb, []))
        tasks.append(tf.transrate_task(assemb,fq1,fq2,unpaired, tf.GEN_PATH_FILTER_FILES(), basename(tf.GEN_PATH_TRANSRATE_DIR()), int(round(float(cpu), 4)), [], transrate_ref)
        for busco_ref in busco_refs:
            tasks.append(tf.busco_task(tf.GEN_PATH_FILTER_FILES(), tf.GEN_PATH_ASSEMBLY(), busco_ref, int(cpu/2), []))
        if(cegma_flag):
             tasks.append(tf.cegma_task(assemb, cpu, []))
    return Supervisor(tasks=tasks,dependencies=dependency_set)

if(__name__ == '__main__'):
    pass
