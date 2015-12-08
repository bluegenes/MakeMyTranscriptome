import argparse
import os
import sys
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf
import time

def gen_quality_supervisor(fasta, transrate_fastq1, transrate_fastq2, transrate_unpaired, busco_ref, cpu=12, cegma_flag=False, transrate_reference='', dependency_set):
    tasks = []
    cegma = tf.cegma_task(cpu, []) 
    busco = tf.busco_task(busco_ref, int(cpu/2), [])
    assembly_stats = tf.assembly_stats_task([])
    transrate = tf.transrate_task(transrate_fastq1, transrate_fastq2, transrate_unpaired, transrate_reference, int(round(float(cpu), 4)), [])
    tasks.append(transrate)
    tasks.append(assembly_stats)
    tasks.append(busco)
    if(cegma_flag):
        tasks.append(cegma)
    return Supervisor(tasks=tasks,dependencies=dependency_set)



if(__name__ == '__main__'):
    pass
