import argparse
import os
import sys
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf



def gen_expression_supervisor(fastq1,fastq2,paired_names,unpaired,unpaired_names,cpu,sample_info,model,dependency_set):
	fasta_to_bed = tf.assembly_to_bed_task([])
	build_bowtie = tf.build_bowtie_task([])
	express_tasks = []
	bowtie_e_tasks = []
	bowtie_i_tasks = []
	sam_sort_tasks = []
	intersect_tasks = []
	for i in range(len(fastq1)):
		bowtie_e = tf.bowtie2_task(fastq1[i],fastq2[i],paired_names[i]+'_express_bt2',0,int(cpu/2),[build_bowtie])
		express = tf.express_task(bowtie_e.targets[0],paired_names[i],[bowtie_e])
		bowtie_i = tf.bowtie2_task(fastq1[i],fastq2[i],paired_names[i]+'_intersect_bt2',1,int(cpu/2),[build_bowtie])
		sam_sort = tf.sam_sort_task(bowtie_i.targets[0],paired_names[i]+'_intersect_bt2_sorted',[bowtie_i])
		intersect_bed = tf.intersect_bed_task(sam_sort.targets[0],fasta_to_bed.targets[0],paired_names[i],[sam_sort,fasta_to_bed])
		bowtie_e_tasks.append(bowtie_e)
		express_tasks.append(express)
		bowtie_i_tasks.append(bowtie_i)
		sam_sort_tasks.append(sam_sort)
		intersect_tasks.append(intersect_bed)
	for i in range(len(unpaired)):
		bowtie_e = tf.bowtie2_unpaired_task(unpaired[i],unpaired_names[i]+'_express_bt2',0,int(cpu/2),[build_bowtie])
		express = tf.express_task(bowtie_e.targets[0],unpaired_names[i],[bowtie_e])
		bowtie_i = tf.bowtie2_unpaired_task(unpaired[i],unpaired_names[i]+'_intersect_bt2',1,int(cpu/2),[build_bowtie])
		sam_sort = tf.sam_sort_task(bowtie_i.targets[0],unpaired_names[i]+'_intersect_bt2_sorted',[bowtie_i])
		intersect_bed = tf.intersect_bed_task(sam_sort.targets[0],fasta_to_bed.targets[0],unpaired_names[i],[sam_sort,fasta_to_bed])
		bowtie_e_tasks.append(bowtie_e)
		express_tasks.append(express)
		bowtie_i_tasks.append(bowtie_i)
		sam_sort_tasks.append(sam_sort)
		intersect_tasks.append(intersect_bed)
	counts_to_table_express = tf.counts_to_table_task([t.targets[0] for t in express_tasks],'express_counts','--eXpress',express_tasks)
	counts_to_table_intersect = tf.counts_to_table_task([t.targets[0] for t in intersect_tasks],'bed_counts','',intersect_tasks)
	deseq2_express = tf.deseq2_task(counts_to_table_express.targets[0],sample_info,'express',model,[counts_to_table_express])
	deseq2_intersect = tf.deseq2_task(counts_to_table_intersect.targets[0],sample_info,'intersect',model,[counts_to_table_intersect])
	all_tasks = [fasta_to_bed,build_bowtie,counts_to_table_intersect,counts_to_table_express,deseq2_intersect,deseq2_express]
	all_tasks.extend(bowtie_e_tasks+bowtie_i_tasks+express_tasks+sam_sort_tasks+intersect_tasks)
	return Supervisor(tasks=all_tasks,dependencies=dependency_set)


if(__name__=='__main__'):
	pass