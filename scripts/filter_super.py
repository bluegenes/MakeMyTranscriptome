import argparse
import os
from os.path import join,exists
import sys
from tasks_v2 import Supervisor, Task
import functions_general as fg
import time


# dependency_set = trinity. but also transrate, for the transrate_good filtering...
def gen_filter_supervisor(main_path_assembly,main_assembly_name,out_dir, transrate_task, dependency_set, tpm_threshold=1):
    tasks = [] 
    filter_full = fg.filter_task(main_path_assembly,main_assembly_name,out_dir,[transrate_task.targets[2]],tpm_threshold,2,[transrate_task])
    tasks.append(filter_full)
    tasks.append(fan.gene_trans_map_task(filter_full.targets[0], out_dir,[filter_full]))
    cp_tr_good = fg.cp_assembly_task(join(out_dir,'good.'+main_assembly_name),transrate_task.targets[1], [transrate_task])
    tasks.append(cp_tr_good)
    tasks.apend(fan.gene_trans_map_task(cp_tr_good.targets[0], out_dir,[cp_tr_good]))
    filter_good = fg.filter_task(cp_tr_good.targets[0],'good.'+ main_assembly_name, out_dir, [transrate_task.targets[2]], tpm_threshold, 2,[cp_tr_good])
    tasks.append(filter_good)
    tasks.append(fan.gene_trans_map_task(filter_good.targets[0], out_dir,[filter_good]))
    return Supervisor(tasks=tasks,dependencies=dependency_set)

if(__name__ == '__main__'):
    pass
