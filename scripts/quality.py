import argparse
import os
from os.path import basename
import sys
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf
import time

def gen_quality_supervisor(transrate_fastq1, transrate_fastq2, transrate_unpaired,  dependency_set, busco_refs, cpu=12, cegma_flag=False, transrate_reference=''):
    tasks = []
    cegma = tf.cegma_task(tf.GEN_PATH_ASSEMBLY(), cpu, []) 
    for busco_ref in busco_refs:
        tasks.append(tf.busco_task(tf.GEN_PATH_QUALITY_FILES(),tf.GEN_PATH_ASSEMBLY(),busco_ref, int(cpu/2), []))
    assembly_stats = tf.assembly_stats_task(tf.GEN_PATH_ASSEMBLY(), [])
    transrate = tf.transrate_task(tf.GEN_PATH_ASSEMBLY(),transrate_fastq1, transrate_fastq2, transrate_unpaired, tf.GEN_PATH_QUALITY_DIR(), basename(tf.GEN_PATH_TRANSRATE_DIR()), int(round(float(cpu), 4)), [], transrate_reference)
    reference_name = os.path.basename(transrate_reference)
    tasks.append(transrate)
    tasks.append(assembly_stats)
    if(cegma_flag):
        tasks.append(cegma)
    return Supervisor(tasks=tasks,dependencies=dependency_set)

if(__name__ == '__main__'):
    pass
