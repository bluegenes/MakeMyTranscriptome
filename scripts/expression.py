import argparse
import os
from os.path import join, dirname, basename
import sys
from tasks_v2 import Supervisor, Task
import functions_general as fg
import functions_annotater as fan
import functions_expression as fex

def gen_expression_supervisor(fastq1,fastq2,paired_names,unpaired,unpaired_names,cpu,sample_info,model,gene_trans_map,dependency_set,assembly_name, assembly_path, out_dir,run_intersectbed=False,run_express=False):
    salmon_tasks, express_tasks, bowtie_e_tasks, bowtie_i_tasks,sam_sort_tasks, intersect_tasks = [],[],[],[],[],[]
    fasta_to_bed = fan.assembly_to_bed_task(assembly_path, out_dir,[])
    build_salmon = fex.build_salmon_task(assembly_path, assembly_name, out_dir,fg.cpumod(cpu, 2),[])
    salmon_gene_map = fex.salmon_gene_map_task(out_dir,assembly_name,gene_trans_map,[])
    if run_express or run_intersectbed:
        build_bowtie = fex.build_bowtie_task(assembly_path,assembly_name, out_dir,[])
        bowtie2_index = join(dirname(build_bowtie.targets[0]), basename(build_bowtie.targets[0]).split('.')[0])
    for i in range(len(fastq1)):
        salmon = fex.salmon_task(build_salmon.targets[0],fastq1[i],fastq2[i],paired_names[i]+'_salmon_'+ assembly_name, salmon_gene_map.targets[0],out_dir,fg.cpumod(cpu,2),[build_salmon, salmon_gene_map])
        salmon_tasks.append(salmon)
        if run_express:
            bowtie_e = fex.bowtie2_task(bowtie2_index,out_dir,fastq1[i],fastq2[i],paired_names[i]+'_express_bt2_'+ assembly_name,0,fg.cpumod(cpu,2),[build_bowtie])
            express = fex.express_task(assembly_path,out_dir,paired_names[i],bowtie_e.targets[0],[bowtie_e])
            bowtie_e_tasks.append(bowtie_e)
            express_tasks.append(express)
        if run_intersectbed:
            bowtie_i = fex.bowtie2_task(bowtie2_index,out_dir,fastq1[i],fastq2[i],paired_names[i]+'_intersect_bt2_'+ assembly_name,1,fg.cpumod(cpu,2),[build_bowtie])
            sam_sort = fex.sam_sort_task(out_dir,bowtie_i.targets[0],paired_names[i]+'_intersect_bt2_sorted_'+ assembly_name,[bowtie_i])
            intersect_bed = fex.intersect_bed_task(out_dir,sam_sort.targets[0],fasta_to_bed.targets[0],paired_names[i],[sam_sort,fasta_to_bed])
            bowtie_i_tasks.append(bowtie_i)
            sam_sort_tasks.append(sam_sort)
            intersect_tasks.append(intersect_bed)
    for i in range(len(unpaired)):
        salmon = fex.salmon_unpaired_task(build_salmon.targets[0],unpaired[i],unpaired_names[i]+'_salmon_'+ assembly_name,salmon_gene_map.targets[0],out_dir,fg.cpumod(cpu,2),[build_salmon,salmon_gene_map])
        salmon_tasks.append(salmon)
        if run_express:
            bowtie_e = fex.bowtie2_unpaired_task(bowtie2_index,out_dir,unpaired[i],unpaired_names[i]+'_express_bt2_'+ assembly_name,0,fg.cpumod(cpu,2),[build_bowtie])
            bowtie_e_tasks.append(bowtie_e)
            express = fex.express_task(assembly_path,out_dir,unpaired_names[i],bowtie_e.targets[0],[bowtie_e])
            express_tasks.append(express)
        if run_intersectbed:
            bowtie_i = fex.bowtie2_unpaired_task(bowtie2_index,out_dir,unpaired[i],unpaired_names[i]+'_intersect_bt2_'+ assembly_name,1,fg.cpumod(cpu,2),[build_bowtie])
            bowtie_i_tasks.append(bowtie_i)
            sam_sort = fex.sam_sort_task(out_dir,bowtie_i.targets[0],unpaired_names[i]+'_intersect_bt2_sorted_'+ assembly_name,[bowtie_i])
            sam_sort_tasks.append(sam_sort)
            intersect_bed = fex.intersect_bed_task(out_dir,sam_sort.targets[0],fasta_to_bed.targets[0],unpaired_names[i],[sam_sort,fasta_to_bed])
            intersect_tasks.append(intersect_bed)
    counts_to_table_salmon=fex.counts_to_table_task(assembly_name, gene_trans_map,out_dir, [t.targets[0] for t in salmon_tasks],assembly_name +'_salmon_counts','--salmon',salmon_tasks)
    deseq2_salmon = fex.deseq2_task(assembly_name,out_dir,counts_to_table_salmon.targets[0],sample_info,assembly_name + '_salmon',model,[counts_to_table_salmon])
    deseq2_salmon_gene = fex.deseq2_task(assembly_name,out_dir,counts_to_table_salmon.targets[1],sample_info,'salmon_gene_'+ assembly_name,model,[counts_to_table_salmon])
    # actually add tasks to supervisor #
    all_tasks = [fasta_to_bed,build_salmon,salmon_gene_map,counts_to_table_salmon, deseq2_salmon, deseq2_salmon_gene]
    all_tasks.extend(salmon_tasks)
    if run_express or run_intersectbed:
        all_tasks.append(build_bowtie)
        if run_express:
            counts_to_table_express = fex.counts_to_table_task(assembly_name,gene_trans_map,out_dir, [t.targets[0] for t in express_tasks],assembly_name +'_express_counts','--eXpress',express_tasks)
            deseq2_express = fex.deseq2_task(assembly_name,out_dir,counts_to_table_express.targets[0],sample_info,'express_'+ assembly_name,model,[counts_to_table_express])
            deseq2_express_gene = fex.deseq2_task(assembly_name,out_dir,counts_to_table_express.targets[1],sample_info,'express_gene_'+ assembly_name,model,[counts_to_table_express])
            all_tasks.extend([counts_to_table_express,deseq2_express, deseq2_express_gene]+bowtie_e_tasks+express_tasks)
        if run_intersectbed:
            counts_to_table_intersect=fex.counts_to_table_task(assembly_name,gene_trans_map,out_dir,[t.targets[0] for t in intersect_tasks],assembly_name +'_bed_counts','',intersect_tasks)
            deseq2_intersect = fex.deseq2_task(assembly_name,out_dir,counts_to_table_intersect.targets[0],sample_info,'intersect_'+ assembly_name,model,[counts_to_table_intersect])
            deseq2_intersect_gene = fex.deseq2_task(assembly_name,out_dir,counts_to_table_intersect.targets[1],sample_info,'intersect_gene_'+ assembly_name,model,[counts_to_table_intersect])
            all_tasks.extend([counts_to_table_intersect,deseq2_intersect, deseq2_intersect_gene]+bowtie_i_tasks+sam_sort_tasks+intersect_tasks)
    return Supervisor(tasks=all_tasks,dependencies=dependency_set)


if(__name__=='__main__'):
        pass
