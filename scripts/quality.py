import argparse
import os
from os.path import basename,join
import sys
from tasks_v2 import Supervisor, Task
import functions_general as fg
import functions_quality as fq
import time

def gen_quality_supervisor(opc, dbs, transrate_fq1, transrate_fq2, dependency_set, busco_refs, assembly_name, assembly_path, out_dir, transrate_dir, reads_dir, filter_dir, cp_transrate=True, cpu=12, cegma_flag=False, transrate_ref=''):
    tasks = []
    for busco_ref in busco_refs:
        tasks.append(fq.busco_task(assembly_path, assembly_name, out_dir, busco_ref, int(cpu/2), []))
    assembly_stats = fq.assembly_stats_task(out_dir,assembly_path, [])
    if transrate_fq1 == None:
        transrate_fq1 = []
    if transrate_fq2 == None:
        transrate_fq2 = []
    transrate = fq.transrate_task(reads_dir,assembly_path,assembly_name,transrate_fq1,transrate_fq2,out_dir,transrate_dir,int(round(float(cpu),4)),[],transrate_ref)
    tasks.append(transrate)
    tasks.append(assembly_stats)
    if cp_transrate:
        tasks.append(fg.cp_assembly_task(join(filter_dir,'good.'+assembly_name),transrate.targets[1], [transrate]))
#    for busco_ref in busco_refs:
#        tasks.append(fq.busco_task(transrate.targets[1], os.path.basename(transrate.targets[1]), out_dir, busco_ref, int(cpu/2), [transrate]))
    if(cegma_flag):
        cegma = fq.cegma_task(out_dir,assembly_path, cpu, []) 
        tasks.append(cegma)
    return Supervisor(tasks=tasks,dependencies=dependency_set)

if(__name__ == '__main__'):
    pass
