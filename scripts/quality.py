import argparse
import os
from os.path import basename
import sys
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf
import time

def gen_quality_supervisor(transrate_fq1, transrate_fq2, transrate_unp, dependency_set, busco_refs, cpu=12, cegma_flag=False, transrate_ref='', assembly_name=tf.NAME_ASSEMBLY, assembly_path= tf.GEN_PATH_ASSEMBLY(), out_dir=tf.GEN_PATH_QUALITY_FILES(), transrate_dir=tf.GEN_PATH_TRANSRATE_DIR(), filter_flag=True, tpm_threshold=0.5):
    tasks = []
    cegma = tf.cegma_task(out_dir,assembly_path, cpu, []) 
    for busco_ref in busco_refs:
        tasks.append(tf.busco_task(assembly_path, assembly_name, out_dir, busco_ref, int(cpu/2), []))
    assembly_stats = tf.assembly_stats_task(out_dir,assembly_path, [])
    transrate = tf.transrate_task(assembly_path,assembly_name,transrate_fq1,transrate_fq2,transrate_unp,out_dir,transrate_dir,int(round(float(cpu), 4)), [], transrate_ref)
    tasks.append(transrate)
    tasks.append(assembly_stats)
    if(cegma_flag):
        tasks.append(cegma)
    transrate_good_assembly = transrate.targets[1]
    if(filter_flag):
        filter_full = tf.filter_task(assembly_path, assembly_name, out_dir, [transrate.targets[2]], tpm_threshold, 2,[transrate])
        filter_good = tf.filter_task(transrate_good_assembly,'good.'+ assembly_name, out_dir, [transrate.targets[2]], tpm_threshold, 2,[transrate])
	tasks.extend([filter_full,filter_good])
        for busco_ref in busco_refs:
            tasks.append(tf.busco_task(filter_full.targets[0],basename(filter_full.targets[0]).split('.fa')[0],out_dir,busco_ref,int(cpu/2),[filter_full]))
            tasks.append(tf.busco_task(filter_good.targets[0],basename(filter_good.targets[0]).split('.fa')[0],out_dir,busco_ref,int(cpu/2),[filter_good]))
    return Supervisor(tasks=tasks,dependencies=dependency_set)

if(__name__ == '__main__'):
    pass
