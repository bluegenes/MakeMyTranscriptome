import argparse
import os
from os.path import basename,join
import sys
from tasks_v2 import Supervisor, Task
import functions_general as fg
import functions_quality as fq
import time

def gen_quality_supervisor(transrate_fq1, transrate_fq2, transrate_unp, dependency_set, busco_refs, cpu=12, cegma_flag=False, transrate_ref='', assembly_name=fg.NAME_ASSEMBLY, assembly_path= fg.GEN_PATH_ASSEMBLY(), out_dir=fg.GEN_PATH_QUALITY_FILES(), transrate_dir=fg.GEN_PATH_TRANSRATE_DIR(), reads_dir=fg.GEN_PATH_ASSEMBLY_FILES()):
    tasks = []
    for busco_ref in busco_refs:
        tasks.append(fq.busco_task(assembly_path, assembly_name, out_dir, busco_ref, int(cpu/2), []))
    assembly_stats = fq.assembly_stats_task(out_dir,assembly_path, [])
    if transrate_fq1 == None:
        transrate_fq1 = []
    if transrate_fq2 == None:
        transrate_fq2 = []
    if transrate_unp == None:
        transrate_unp = []
    transrate = fq.transrate_task(reads_dir,assembly_path,assembly_name,transrate_fq1,transrate_fq2,transrate_unp,out_dir,transrate_dir,int(round(float(cpu),4)),[],transrate_ref)
    tasks.append(transrate)
    tasks.append(assembly_stats)
    busco_transrate_good = fq.busco_task(transrate.targets[1], os.path.basename(transrate.targets[1]), out_dir, busco_ref, int(cpu/2), [transrate])
    tasks.append(busco_transrate_good)
    if(cegma_flag):
        cegma = fq.cegma_task(out_dir,assembly_path, cpu, []) 
        tasks.append(cegma)
    return Supervisor(tasks=tasks,dependencies=dependency_set)

if(__name__ == '__main__'):
    pass
