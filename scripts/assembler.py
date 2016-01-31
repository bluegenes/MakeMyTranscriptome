import argparse
import os
import sys
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf
import time


def gen_unpaired_prinseq_supervisor(fastq1, fastq2, unpaired, dependency_set, rmdup):
    tasks = []
    prinseq_count = len(fastq1)
    prinseq_opts = '--derep 14' if(rmdup) else ''
    for input1 in unpaired:
        p_task = tf.prinseq_unpaired_task(input1,'prinseq_output_'+str(prinseq_count),prinseq_opts,[])
        prinseq_count += 1
        tasks.append(p_task)
    return Supervisor(tasks=tasks)


def gen_paired_prinseq_supervisor(fastq1,fastq2,unpaired,dependency_set,rmdup):
    tasks = []
    prinseq_count = 0
    prinseq_opts = '--derep 14' if(rmdup) else ''
    for input1,input2 in zip(fastq1,fastq2):
        p_task = tf.prinseq_task(input1, input2, 'prinseq_output_'+str(prinseq_count), prinseq_opts, [])
        prinseq_count += 1
        tasks.append(p_task)
    return Supervisor(tasks=tasks)


def gen_unpaired_trimmomatic_supervisor(fq1, fq2, unpaired, dependency_set, cpu_cap):
    tasks = []
    count = len(fq1)
    cpu_mod = int(round(float(cpu_cap)/len(unpaired)))
    for i in unpaired:
        trim_task = tf.trimmomatic_unpaired_task(i, cpu_mod, 'trimmomatic_output_'+str(count), dependency_set)
        count+=1
        tasks.append(trim_task)
    return Supervisor(tasks=tasks)        


def gen_paired_trimmomatic_supervisor(fq1, fq2, unpaired, dependency_set, cpu_cap):
    tasks = []
    count = 0
    cpu_mod = int(round(float(cpu_cap)/len(fq1)))
    for i1, i2 in zip(fq1, fq2):
        trim_task = tf.trimmomatic_task(i1, i2, cpu_mod, 'trimmomatic_output_'+str(count), dependency_set)
        count += 1
        tasks.append(trim_task)
    return Supervisor(tasks=tasks)


#def gen_assembly_supervisor(fastq1, fastq2, unpaired, dependency_set, busco_ref, no_trim=False, rnaSPAdes=False, rmdup=False, subset_size=50000000, cpu=12, cegma_flag=False, subset_seed='I am a seed value', normalize_flag=False, truncate_opt=-1, transrate_reference='', trimmomatic_flag=False):
def gen_assembly_supervisor(fastq1, fastq2, unpaired, dependency_set, no_trim=False, rnaSPAdes=False, rmdup=False, subset_size=50000000, cpu=12, subset_seed='I am a seed value', normalize_flag=False, truncate_opt=-1, trimmomatic_flag=False):
    tasks = []
    tasks.append(tf.fastqc_task(fastq1+fastq2+unpaired, 'pre_trimming', []))
    assembler_dependencies = []
    transrate_fastq1 = fastq1
    transrate_fastq2 = fastq2
    transrate_unpaired = unpaired
    if(fastq1 != []):
        if(not no_trim):
            if(trimmomatic_flag):
                paired_sup = gen_paired_trimmomatic_supervisor(fastq1, fastq2, unpaired, [], cpu)
            else:
                paired_sup = gen_paired_prinseq_supervisor(fastq1, fastq2, unpaired, [], rmdup)
            fastq1 = [paired_sup.targets[x] for x in range(0, len(paired_sup.targets), 2)]
            fastq2 = [paired_sup.targets[x] for x in range(1, len(paired_sup.targets), 2)]
            transrate_fastq1 = fastq1
            transrate_fastq2 = fastq2
            tasks.append(paired_sup)
            tasks.append(tf.fastqc_task(fastq1+fastq2, 'post_trimming_paired', [paired_sup]))
        subset_dependencies = [paired_sup] if(not no_trim) else []
        subset = tf.subset_task(fastq1, fastq2, 'final_reads', subset_size, subset_seed, subset_dependencies)
        fastq1 = [subset.targets[0]]
        fastq2 = [subset.targets[1]]
        tasks.append(subset)
        assembler_dependencies = [subset]
        late_fastqc = tf.fastqc_task(subset.targets, 'final_reads_paired', [subset])
        tasks.append(late_fastqc)
        if(truncate_opt >= 0):
            truncate = tf.truncate_task(fastq1[0], fastq2[0], truncate_opt, [subset])
            fastq1 = [truncate.targets[0]]
            fastq2 = [truncate.targets[1]]
            assembler_dependencies = [truncate]
            tasks.append(truncate)
    if(unpaired != []):
        if(not no_trim):
            if(trimmomatic_flag):
                unpaired_sup = gen_unpaired_prinseq_supervisor(fastq1, fastq2, unpaired, [], cpu)
            else:
                unpaired_sup = gen_unpaired_prinseq_supervisor(fastq1, fastq2, unpaired, [], rmdup)
            unpaired = unpaired_sup.targets
            transrate_unpaired = unpaired
            tasks.append(unpaired_sup)
            tasks.append(tf.fastqc_task(unpaired, 'post_trimming_unpaired', [unpaired_sup]))
            assembler_dependencies.append(unpaired_sup)
    if(rnaSPAdes):
        rnaspades = tf.rnaspades_task(fastq1, fastq2, unpaired, cpu, assembler_dependencies)
        tasks.append(rnaspades)
    else:
        trinity = tf.trinity_task(fastq1, fastq2, unpaired, cpu, int(cpu/2), 120, 120, normalize_flag, assembler_dependencies)
        tasks.append(trinity)
    assembler_main_task = tasks[-1]
    return Supervisor(tasks=tasks)


def run_assembler():
    pass


if(__name__ == '__main__'):
    pass
