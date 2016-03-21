import argparse
import os
import sys
from tasks_v2 import Supervisor, Task
import functions_general as fg
import functions_assembler as fa
import functions_annotater as fan
import time


def gen_unpaired_prinseq_supervisor(out_dir,fastq1, fastq2, unpaired, dependency_set, rmdup):
    tasks = []
    prinseq_count = len(fastq1)
    prinseq_opts = '--derep 14' if(rmdup) else ''
    for input1 in unpaired:
        p_task = fa.prinseq_unpaired_task(out_dir,input1,'prinseq_output_'+str(prinseq_count),prinseq_opts,[])
        prinseq_count += 1
        tasks.append(p_task)
    return Supervisor(tasks=tasks)


def gen_paired_prinseq_supervisor(out_dir,fastq1,fastq2,unpaired,dependency_set,rmdup):
    tasks = []
    prinseq_count = 0
    prinseq_opts = '--derep 14' if(rmdup) else ''
    for input1,input2 in zip(fastq1,fastq2):
        p_task = fa.prinseq_task(out_dir,input1, input2, 'prinseq_output_'+str(prinseq_count), prinseq_opts, [])
        prinseq_count += 1
        tasks.append(p_task)
    return Supervisor(tasks=tasks)


def gen_unpaired_trimmomatic_supervisor(out_dir,fq1, fq2, unpaired, dependency_set, cpu_cap):
    tasks = []
    count = len(fq1)
#    cpu_mod = min(len(fq1),cpu_cap)
    cpu_mod = int(round(float(cpu_cap)/len(unpaired)))
    for i in unpaired:
        trim_task = fa.trimmomatic_unpaired_task(out_dir,i, cpu_mod, 'trimmomatic_output_'+str(count), dependency_set)
        count+=1
        tasks.append(trim_task)
    return Supervisor(tasks=tasks)        


def gen_paired_trimmomatic_supervisor(out_dir,fq1, fq2, unpaired, dependency_set, cpu_cap):
    tasks = []
    count = 0
#    cpu_mod = min(len(fq1),cpu_cap) 
    cpu_mod = int(round(float(cpu_cap)/len(fq1)))
    for i1, i2 in zip(fq1, fq2):
        trim_task = fa.trimmomatic_task(out_dir,i1, i2, cpu_mod, 'trimmomatic_output_'+str(count), dependency_set)
        count += 1
        tasks.append(trim_task)
    return Supervisor(tasks=tasks)

def gen_trimming_supervisor(out_dir,fq1,fq2,unpaired,no_trim,trimmomatic_flag,rmdup,subset_size,subset_seed,truncate_opt,dependency_set,cpu_cap):
    tasks = []
    deps = []
    if (not no_trim):
        if(fq1 != []):
            if(trimmomatic_flag):
                paired_sup = gen_paired_trimmomatic_supervisor(out_dir,fq1, fq2, unpaired, dependency_set, cpu_cap)
            else:
                paired_sup = gen_paired_prinseq_supervisor(out_dir,fq1, fq2, unpaired,  dependency_set, rmdup)
            fq1 = [paired_sup.targets[x] for x in range(0, len(paired_sup.targets), 2)]
            fq2 = [paired_sup.targets[x] for x in range(1, len(paired_sup.targets), 2)]
            tasks.append(paired_sup)
            tasks.append(fa.fastqc_task(out_dir, fq1+fq2, 'post_trimming_paired',int(round(float(cpu_cap)/2)),[paired_sup]))
	    deps.append(paired_sup)
        if(unpaired != []):
            if(trimmomatic_flag):
                unpaired_sup = gen_unpaired_trimmomatic_supervisor(out_dir,fq1,fq2, unpaired, dependency_set, cpu_cap)
            else:
                unpaired_sup = gen_unpaired_prinseq_supervisor(out_dir,fq1,fq2,unpaired,dependency_set,rmdup)
            unpaired = unpaired_sup.targets
            tasks.append(unpaired_sup)
            tasks.append(fa.fastqc_task(out_dir,unpaired, 'post_trimming_unpaired',int(round(float(cpu_cap)/2)),[unpaired_sup]))
            deps.append(unpaired_sup)
    subset = fa.subset_task(out_dir, fq1, fq2, 'final_reads', subset_size, subset_seed, deps)
    fq1 = [subset.targets[0]]
    fq2 = [subset.targets[1]]
    tasks.append(subset)
    if(truncate_opt >= 0):
        truncate = fa.truncate_task(out_dir, fastq1[0], fastq2[0], truncate_opt, [subset])
        fq1 = [truncate.targets[0]]
        fq2 = [truncate.targets[1]]
        deps.append(truncate)
        tasks.append(truncate)
    if any([truncate_opt>=0, subset_size < 10**15]): # some subsetting may have occurred
        late_fastqc = fa.fastqc_task(out_dir, subset.targets, 'final_reads_paired',cpu_cap,deps)
        tasks.append(late_fastqc)
    return Supervisor(tasks=tasks)

def gen_assembly_supervisor(out_dir, fastq1, fastq2, unpaired, dependency_set, no_trim=False, rnaSPAdes=False, rmdup=False, subset_size=50000000, cpu=12, subset_seed='I am a seed value', normalize_flag=False, truncate_opt=-1, trimmomatic_flag=True, path_assembly=fg.GEN_PATH_ASSEMBLY()):
    trinity_memory = 160 # make this a user option
    tasks = []
    tasks.append(fa.fastqc_task(out_dir,fastq1+fastq2+unpaired,'pre_trimming',min(cpu,len(fastq1+fastq2+unpaired)), []))
    trim_reads = gen_trimming_supervisor(out_dir,fastq1,fastq2,unpaired,no_trim,trimmomatic_flag,rmdup,subset_size,subset_seed, truncate_opt,[],cpu)
    tasks.append(trim_reads)
    if(rnaSPAdes):
        rnaspades = fa.rnaspades_task(path_assembly, out_dir,fastq1, fastq2, unpaired, cpu, [trim_reads])
        tasks.append(rnaspades)
    else:
        trinity = fa.trinity_task(path_assembly, out_dir, fastq1, fastq2, unpaired, cpu, int(cpu/2), trinity_memory, trinity_memory, normalize_flag, [trim_reads])
        tasks.append(trinity)
        gene_trans_map = fan.gene_trans_map_task(path_assembly,out_dir,[trinity])
        tasks.append(gene_trans_map)
    return Supervisor(tasks=tasks)


def run_assembler():
    pass


if(__name__ == '__main__'):
    pass
