import argparse
import os
from os.path import join,exists
import sys
from tasks_v2 import Supervisor, Task
import functions_general as fg
import time


def gen_filter_supervisor(main_path_assembly,main_assembly_name,out_dir, transrate_task, dependency_set, tpm_threshold=1):
    tasks = [] 
    filter_full = fg.filter_task(main_path_assembly,main_assembly_name,out_dir,[transrate_task.targets[2]],tpm_threshold,2,[transrate_task])
    tasks.append(filter_full)
    filter_good = fg.filter_task(transrate_task.targets[1],'good.'+ main_assembly_name, out_dir, [transrate_task.targets[2]], tpm_threshold, 2,[transrate_task])
    tasks.append(filter_good)
    return Supervisor(tasks=tasks,dependencies=dependency_set)

if(__name__ == '__main__'):
    pass
