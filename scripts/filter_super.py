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
    cp_tr_good = fg.cp_assembly_task(join(out_dir,'good.'+main_assembly_name),transrate_task.targets[1], [transrate_task])
    tasks.append(cp_tr_good)
    filter_good = fg.filter_task(cp_tr_good.targets[0],'good.'+ main_assembly_name, out_dir, [transrate_task.targets[2]], tpm_threshold, 2,[cp_tr_good])
    tasks.append(filter_good)


    #filter full assembly by tpm; assess
    #filter_full = fg.filter_task(assembly_path,assembly_name,out_dir,[transrate_task.targets[2]],tpm_threshold,2,[transrate_task])
    #tasks.append(filter_full)
   # tasks.append(fg.assembly_stats_task(out_dir,filter_full.targets[0], [filter_full]))
   # tasks.append(fg.transrate_task(filter_full.targets[0],fq1,fq2,unpaired,out_dir,basename(transrate_dir),int(round(float(cpu), 4)),[filter_full],transrate_ref)
    

    #for busco_ref in busco_refs:
    #    tasks.append(fg.busco_task(filter_full.targets[0],basename(filter_full.targets[0]).split('.fa')[0],out_dir,busco_ref,int(cpu/2),[filter_full]))
    #expression: filtered_full
    #filter_full_build_salmon = fg.build_salmon_task(filter_full.targets[0],basename(filter_full.targets[0]).split('.fa')[0],out_dir,cpu,[filter_full])
    #assess the transrate_good assembly:
    #cp_tr_good = fg.cp_assembly_task(join(filter_files,'good.'+assembly_name),transrate_task.targets[1], [transrate_task])
    #tasks.append(cp_tr_good)
    #tasks.append(fg.assembly_stats_task(out_dir,cp_tr_good.targets[0], [cp_tr_good]))
    

    #for busco_ref in busco_refs:
    #    tasks.append(fg.busco_task(cp_tr_good.targets[0],basename(cp_tr_good.targets[0]).split('.fa')[0],out_dir,busco_ref,int(cpu/2),[cp_tr_good]))
    #filter transrate good assembly: assess
    #filter_good = fg.filter_task(cp_tr_good.targets[0],'good.'+ assembly_name, out_dir, [transrate_task.targets[2]], tpm_threshold, 2,[cp_tr_good])
    #tasks.append(filter_good)
    #tasks.append(fg.assembly_stats_task(out_dir,filter_good.targets[0], [filter_good]))
    #tasks.append(fg.transrate_task(filter_good.targets[0],fq1,fq2,unpaired,out_dir,basename(transrate_dir),int(round(float(cpu), 4)),[filter_good],transrate_ref)
    #for busco_ref in busco_refs:
    #    tasks.append(fg.busco_task(filter_good.targets[0],basename(filter_good.targets[0]).split('.fa')[0],out_dir,busco_ref,int(cpu/2),[filter_good]))
    #salmon_transrate_good = fg.salmon_task()
    #salmon_full = fg.salmon_task()
    return Supervisor(tasks=tasks,dependencies=dependency_set)

if(__name__ == '__main__'):
    pass
