import argparse
import os
from os.path import join, dirname, basename
import sys
from tasks_v2 import Supervisor, Task
import functions_general as fg
import functions_annotater as fan
import functions_expression as fex
import assembler as assemb

salmon_naming = 'salmon'
express_naming = 'express'
intersect_naming = 'intersect'

def gen_salmon_supervisor(opc, dbs, fastq1,fastq2,paired_names,unpaired,unpaired_names,assembly_path,assembly_name,gene_trans_map,sample_info,model,out_dir,cpu_cap, deps):
    salmon_tasks = []
    salmon_dir = fg.make_dir_task(os.path.join(out_dir,'salmon'))
    out_dir = salmon_dir.targets[0]
    build_salmon = fex.build_salmon_task(assembly_path, assembly_name, out_dir,fg.cpumod(cpu_cap, 2),[salmon_dir])
    salmon_gene_map = fex.salmon_gene_map_task(out_dir,assembly_name,gene_trans_map,[salmon_dir])
    deps = deps + [build_salmon, salmon_gene_map]
    for i in range(len(fastq1)):
        #filename = '_'.join([paired_names[i],salmon_naming,assembly_name]) 
        filename = paired_names[i] #,salmon_naming,assembly_name]) 
        salmon = fex.salmon_task(build_salmon.targets[0],fastq1[i],fastq2[i],filename, salmon_gene_map.targets[0],out_dir,fg.cpumod(cpu_cap,2),deps)
        salmon_tasks.append(salmon)
    for i in range(len(unpaired)):
        #filename = '_'.join([unpaired_names[i],salmon_naming,assembly_name]) 
        filename = unpaired_names[i] #,salmon_naming,assembly_name]) 
        salmon = fex.salmon_unpaired_task(build_salmon.targets[0],unpaired[i],filename,salmon_gene_map.targets[0],out_dir,fg.cpumod(cpu_cap,2),deps)
        salmon_tasks.append(salmon)
    transcriptName = assembly_name  #'_'.join([assembly_name,salmon_naming])
    geneName = assembly_name + '_gene' #'_'.join([assembly_name,salmon_naming,'gene'])
    counts_to_table_salmon=fex.counts_to_table_task(assembly_name,gene_trans_map,out_dir,[t.targets[0] for t in salmon_tasks],transcriptName,'--salmon',salmon_tasks)
    deseq2_salmon = fex.deseq2_task(assembly_name,out_dir,counts_to_table_salmon.targets[0],sample_info,transcriptName,model,[counts_to_table_salmon])
    deseq2_salmon_gene = fex.deseq2_task(assembly_name,out_dir,counts_to_table_salmon.targets[1],sample_info,geneName,model,[counts_to_table_salmon])
    salmon_tasks = [salmon_dir,build_salmon,salmon_gene_map,counts_to_table_salmon, deseq2_salmon, deseq2_salmon_gene]+salmon_tasks
    return Supervisor(tasks = salmon_tasks)

def gen_express_supervisor(fastq1,fastq2,paired_names,unpaired,unpaired_names,assembly_path,assembly_name,bowtie2_index,gene_trans_map,sample_info,model,out_dir,cpu_cap,deps):
    express_tasks,bowtie_e_tasks = [],[]
    express_dir = fg.make_dir_task(os.path.join(out_dir,'express'))
    out_dir = express_dir.targets[0]
    for i in range(len(fastq1)):
        filename = paired_names[i] #'_'.join([paired_names[i],express_naming,assembly_name]) 
        #filename = '_'.join([paired_names[i],express_naming,assembly_name]) 
	bowtie_e = fex.bowtie2_task(bowtie2_index,out_dir,fastq1[i],fastq2[i],filename,0,fg.cpumod(cpu_cap,2),deps)
        express = fex.express_task(bowtie2_index,assembly_path,out_dir,paired_names[i],bowtie_e.targets[0],[bowtie_e])
        bowtie_e_tasks.append(bowtie_e)
        express_tasks.append(express)
    for i in range(len(unpaired)):
        filename = unpaired_names[i] #'_'.join([unpaired_names[i],express_naming,assembly_name])
        bowtie_e = fex.bowtie2_unpaired_task(bowtie2_index,out_dir,unpaired[i],filename,0,fg.cpumod(cpu_cap,2),deps)
        bowtie_e_tasks.append(bowtie_e)
        express = fex.express_task(bowtie2_index,assembly_path,out_dir,unpaired_names[i],bowtie_e.targets[0],[bowtie_e])
        express_tasks.append(express)
    transcriptName = assembly_name #'_'.join([assembly_name,express_naming])
    geneName = assembly_name + '_gene' #'_'.join([assembly_name,express_naming,'gene'])
    counts_to_table_express = fex.counts_to_table_task(assembly_name,gene_trans_map,out_dir,[t.targets[0] for t in express_tasks],transcriptName,'--eXpress',express_tasks)
    deseq2_express = fex.deseq2_task(assembly_name,out_dir,counts_to_table_express.targets[0],sample_info,transcriptName,model,[counts_to_table_express])
    deseq2_express_gene = fex.deseq2_task(assembly_name,out_dir,counts_to_table_express.targets[1],sample_info,geneName,model,[counts_to_table_express])
    e_tasks = [express_dir,counts_to_table_express,deseq2_express,deseq2_express_gene]+bowtie_e_tasks+express_tasks
    return Supervisor(tasks = e_tasks)

def gen_intersect_supervisor(fq1,fq2,paired_names,unpaired,unpaired_names,assembly_path,assembly_name,bowtie2_index,gene_trans_map,sample_info,model,out_dir,cpu_cap, deps):
    intersect_tasks,bowtie_i_tasks,sam_sort_tasks = [],[],[]
    intersect_dir = fg.make_dir_task(os.path.join(out_dir,'intersectBed'))
    out_dir = intersect_dir.targets[0]
    deps.append(intersect_dir)
    fasta_to_bed = fan.assembly_to_bed_task(assembly_path, out_dir,[intersect_dir])
    for i in range(len(fq1)):
        filename = paired_names[i] #'_'.join([paired_names[i],intersect_naming,assembly_name]) 
        #filename = '_'.join([paired_names[i],intersect_naming,assembly_name]) 
        bowtie_i = fex.bowtie2_task(bowtie2_index,out_dir,fq1[i],fq2[i],filename,1,fg.cpumod(cpu_cap,2),deps)
	sorted_name = filename + '_sorted'
        sam_sort = fex.sam_sort_task(out_dir,bowtie_i.targets[0],sorted_name,[bowtie_i])
        intersect_bed = fex.intersect_bed_task(out_dir,sam_sort.targets[0],fasta_to_bed.targets[0],paired_names[i],[sam_sort,fasta_to_bed])
        bowtie_i_tasks.append(bowtie_i)
        sam_sort_tasks.append(sam_sort)
        intersect_tasks.append(intersect_bed)
    for i in range(len(unpaired)):
        filename = unpaired_names[i] #'_'.join([unpaired_names[i],intersect_naming,assembly_name])
        bowtie_i = fex.bowtie2_unpaired_task(bowtie2_index,out_dir,unpaired[i],filename,1,fg.cpumod(cpu_cap,2),deps)
        bowtie_i_tasks.append(bowtie_i)
	sorted_name = filename + '_sorted'
        sam_sort = fex.sam_sort_task(out_dir,bowtie_i.targets[0],sorted_name,[bowtie_i])
        sam_sort_tasks.append(sam_sort)
        intersect_bed = fex.intersect_bed_task(out_dir,sam_sort.targets[0],fasta_to_bed.targets[0],unpaired_names[i],[sam_sort,fasta_to_bed])
        intersect_tasks.append(intersect_bed)
    transcriptName = assembly_name #'_'.join([assembly_name,express_naming])
    geneName = assembly_name + '_gene' #'_'.join([assembly_name,express_naming,'gene'])
    counts_to_table_intersect=fex.counts_to_table_task(assembly_name,gene_trans_map,out_dir,[t.targets[0] for t in intersect_tasks],transcriptName,'',intersect_tasks)
    deseq2_intersect = fex.deseq2_task(assembly_name,out_dir,counts_to_table_intersect.targets[0],sample_info,transcriptName,model,[counts_to_table_intersect])
    deseq2_intersect_gene = fex.deseq2_task(assembly_name,out_dir,counts_to_table_intersect.targets[1],sample_info,geneName,model,[counts_to_table_intersect])
    i_tasks = [intersect_dir,fasta_to_bed,counts_to_table_intersect,deseq2_intersect, deseq2_intersect_gene]+bowtie_i_tasks+sam_sort_tasks+intersect_tasks
    return Supervisor(tasks=i_tasks)

def gen_expression_supervisor(fastq1,fastq2,paired_names,unpaired,unpaired_names,cpu,sample_info,model,gene_trans_map,dependency_set,assembly_name, assembly_path, out_dir,run_express=False,run_intersectbed=False):
    all_tasks = []
    deps = []
    trim_reads = False
    if trim_reads:
        trimmomatic_flag = True
	rmdup = False
	truncate_opt = False
	trim_tasks,fastq1,fastq2,unpaired=assemb.gen_trimming_supervisor(out_dir,fastq1,fastq2,unpaired,False,trimmomatic_flag,rmdup,10**15,0,truncate_opt,[],cpu) 
	all_tasks.append(trim_tasks)
	deps.append(trim_tasks)
    salmon_tasks = gen_salmon_supervisor(fastq1,fastq2,paired_names,unpaired,unpaired_names,assembly_path,assembly_name,gene_trans_map,sample_info,model,out_dir,cpu, deps)
    all_tasks.append(salmon_tasks)
    if run_express or run_intersectbed:
        build_bowtie = fex.build_bowtie_task(assembly_path,assembly_name, out_dir,[])
        bowtie2_index = join(dirname(build_bowtie.targets[0]),basename(build_bowtie.targets[0]).split('.')[0])
	all_tasks.append(build_bowtie)
	if run_express:
	    express_tasks = gen_express_supervisor(fastq1,fastq2,paired_names,unpaired,unpaired_names,assembly_path,assembly_name,bowtie2_index,gene_trans_map,sample_info,model,out_dir,cpu, [build_bowtie])
            all_tasks.append(express_tasks)
	if run_intersectbed:
	    intersect_tasks = gen_intersect_supervisor(fastq1,fastq2,paired_names,unpaired,unpaired_names,assembly_path,assembly_name,bowtie2_index,gene_trans_map,sample_info,model,out_dir,cpu,[build_bowtie])
	    all_tasks.append(intersect_tasks)
    return Supervisor(tasks=all_tasks,dependencies=dependency_set)


if(__name__=='__main__'):
        pass
