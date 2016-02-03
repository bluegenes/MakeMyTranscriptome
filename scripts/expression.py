import argparse
import os
from os.path import join, dirname, basename
import sys
from tasks_v2 import Supervisor, Task
#import task_functions_v2 as tf
import functions_general as fg
import functions_annotater as fan
import functions_expression as fex

def gen_expression_supervisor(fastq1,fastq2,paired_names,unpaired,unpaired_names,cpu,sample_info,model,dependency_set,run_intersectbed=False,run_express=False, assembly_name=fg.NAME_ASSEMBLY, assembly_path= fg.GEN_PATH_ASSEMBLY(), out_dir=fg.GEN_PATH_EXPRESSION_FILES(), gene_trans_map=''):
    fasta_to_bed = fan.assembly_to_bed_task(assembly_path, out_dir,[])
    build_salmon = fex.build_salmon_task(assembly_path, assembly_name, out_dir,cpu,[])
#    gene_trans_map = assembly_path.rsplit('.fa')[0] + '.gene_trans_map'
    salmon_gene_map = fex.salmon_gene_map_task(out_dir,assembly_name,gene_trans_map,[])
    build_bowtie = fex.build_bowtie_task(assembly_path,assembly_name, out_dir,[])
    bowtie2_index = join(dirname(build_bowtie.targets[0]), basename(build_bowtie.targets[0]).split('.')[0])
    salmon_tasks, express_tasks, bowtie_e_tasks, bowtie_i_tasks,sam_sort_tasks, intersect_tasks = [],[],[],[],[],[]
    for i in range(len(fastq1)):
        bowtie_e = fex.bowtie2_task(bowtie2_index,out_dir,fastq1[i],fastq2[i],paired_names[i]+'_express_bt2',0,int(cpu/2),[build_bowtie])
        express = fex.express_task(assembly_path,out_dir,paired_names[i],bowtie_e.targets[0],[bowtie_e])
        bowtie_i = fex.bowtie2_task(bowtie2_index,out_dir,fastq1[i],fastq2[i],paired_names[i]+'_intersect_bt2',1,int(cpu/2),[build_bowtie])
        sam_sort = fex.sam_sort_task(out_dir,bowtie_i.targets[0],paired_names[i]+'_intersect_bt2_sorted',[bowtie_i])
        intersect_bed = fex.intersect_bed_task(out_dir,sam_sort.targets[0],fasta_to_bed.targets[0],paired_names[i],[sam_sort,fasta_to_bed])
        salmon = fex.salmon_task(build_salmon.targets[0],fastq1[i],fastq2[i],paired_names[i]+'_salmon', salmon_gene_map.targets[0],out_dir,int(cpu/2),[build_salmon, salmon_gene_map])
        salmon_tasks.append(salmon)
        if run_express:
            bowtie_e_tasks.append(bowtie_e)
            express_tasks.append(express)
        if run_intersectbed:
            bowtie_i_tasks.append(bowtie_i)
            sam_sort_tasks.append(sam_sort)
            intersect_tasks.append(intersect_bed)
    for i in range(len(unpaired)):
        bowtie_e = fex.bowtie2_unpaired_task(bowtie2_index,out_dir,unpaired[i],unpaired_names[i]+'_express_bt2',0,int(cpu/2),[build_bowtie])
        express = fex.express_task(assembly_path,out_dir,unpaired_names[i],bowtie_e.targets[0],[bowtie_e])
        bowtie_i = fex.bowtie2_unpaired_task(bowtie2_index,out_dir,unpaired[i],unpaired_names[i]+'_intersect_bt2',1,int(cpu/2),[build_bowtie])
        sam_sort = fex.sam_sort_task(out_dir,bowtie_i.targets[0],unpaired_names[i]+'_intersect_bt2_sorted',[bowtie_i])
        intersect_bed = fex.intersect_bed_task(out_dir,sam_sort.targets[0],fasta_to_bed.targets[0],unpaired_names[i],[sam_sort,fasta_to_bed])
        salmon = fex.salmon_unpaired_task(build_salmon.targets[0],unpaired_names[i],unpaired_names[i]+'_salmon',salmon_gene_map.targets[0],out_dir,int(cpu/2),[build_salmon,salmon_gene_map])
        salmon_tasks.append(salmon)
        if run_express:
            bowtie_e_tasks.append(bowtie_e)
            express_tasks.append(express)
        if run_intersectbed:
            bowtie_i_tasks.append(bowtie_i)
            sam_sort_tasks.append(sam_sort)
            intersect_tasks.append(intersect_bed)
    counts_to_table_express = fex.counts_to_table_task(gene_trans_map,out_dir, [t.targets[0] for t in express_tasks],'express_counts','--eXpress',express_tasks)
    counts_to_table_intersect=fex.counts_to_table_task(gene_trans_map,out_dir,[t.targets[0] for t in intersect_tasks],'bed_counts','',intersect_tasks)
    counts_to_table_salmon=fex.counts_to_table_task(gene_trans_map,out_dir, [t.targets[0] for t in salmon_tasks],'salmon_counts','--salmon',salmon_tasks)
    deseq2_express = fex.deseq2_task(out_dir,counts_to_table_express.targets[0],sample_info,'express',model,[counts_to_table_express])
    deseq2_intersect = fex.deseq2_task(out_dir,counts_to_table_intersect.targets[0],sample_info,'intersect',model,[counts_to_table_intersect])
    deseq2_salmon = fex.deseq2_task(out_dir,counts_to_table_salmon.targets[0],sample_info,'salmon',model,[counts_to_table_salmon])

    # actually add tasks to supervisor #
    all_tasks = [fasta_to_bed,build_salmon,gene_trans_map,salmon_gene_map,counts_to_table_salmon, deseq2_salmon]
    all_tasks.extend(salmon_tasks)
    if run_express or run_intersectbed:
        all_tasks.append(build_bowtie)
        if run_express:
            all_tasks.extend([counts_to_table_express,deseq2_express]+bowtie_e_tasks+express_tasks)
        if run_intersectbed:
	    all_tasks.extend([counts_to_table_intersect,deseq2_intersect]+bowtie_i_tasks+sam_sort_tasks+intersect_tasks)
    return Supervisor(tasks=all_tasks,dependencies=dependency_set)


if(__name__=='__main__'):
        pass
