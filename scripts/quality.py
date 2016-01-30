import argparse
import os
from os.path import basename
import sys
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf
import time

def gen_quality_supervisor(transrate_fq1, transrate_fq2, transrate_unp, dependency_set, busco_refs, cpu=12, cegma_flag=False, transrate_ref='', assembly_name=tf.NAME_ASSEMBLY, assembly_path= tf.GEN_PATH_ASSEMBLY(), out_dir=tf.GEN_PATH_QUALITY_FILES(), transrate_dir=tf.GEN_PATH_TRANSRATE_DIR()):
    tasks = []
    cegma = tf.cegma_task(out_dir,assembly_path, cpu, []) 
    for busco_ref in busco_refs:
        tasks.append(tf.busco_task(out_dir, assembly_path, busco_ref, int(cpu/2), []))
    assembly_stats = tf.assembly_stats_task(assembly_path, [])
    transrate = tf.transrate_task(assembly_path,transrate_fq1,transrate_fq2,transrate_unp,out_dir,basename(transrate_dir),int(round(float(cpu), 4)), [], transrate_ref)
    tasks.append(transrate)
    tasks.append(assembly_stats)
    if(cegma_flag):
        tasks.append(cegma)
    return Supervisor(tasks=tasks,dependencies=dependency_set)

if(__name__ == '__main__'):
    pass
