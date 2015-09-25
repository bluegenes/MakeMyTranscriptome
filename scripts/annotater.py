import argparse
import os
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf


def gen_annotation_assembly(cpu,busco_ref,blast_uniref90,cegma,dependency_set):
	cegma = tf.cegma_task(cpu,[])
	busco = tf.busco_task(busco_ref,int(cpu/2),[])
	assembly_stats = tf.assembly_stats_task([])
	gene_trans_map = tf.gene_trans_map_task([])
	blastx_sprot = tf.blastx_task(tf.PATH_SWISS_PROT,int(cpu/2),[])
	rnammer = tf.rnammer_task([])
	predict_orfs = tf.predict_orfs_task(int(round(cpu/2)),[])
	pep_path = predict_orfs.targets[0]
	signalp = tf.signalp_task(pep_path,[predict_orfs])
	blastp_sprot = tf.blastp_task(pep_path,tf.PATH_SWISS_PROT,int(cpu/2),[predict_orfs])
	tmhmm = tf.tmhmm_task(pep_path,[predict_orfs])
	pfam = tf.pfam_task(pep_path,int(cpu/2),[predict_orfs])
	dependencies = [gene_trans_map,blastx_sprot,rnammer,predict_orfs,
					blastp_sprot,pfam,signalp,tmhmm]
	if(blast_uniref90):
		blastp_ur90 = tf.blastp_task(pep_path,tf.PATH_UNIREF90,int(cpu/2),[predict_orfs])
		blastx_ur90 = tf.blastx_task(tf.PATH_UNIREF90,int(cpu)/2,[])
		dependencies.append(blastx_ur90)
		dependencies.append(blastp_ur90)
		blastp_ur90_target = blastp_ur90.targets[0]
		blastx_ur90_target = blastx_ur90.targets[0]
	else:
		blastp_ur90_target = 'NONE'
		blastx_ur90_target = 'NONE'
	annot = tf.annot_table_task(gene_trans_map.targets[0],blastx_sprot.targets[0],
		blastx_ur90_target, rnammer.targets[0], predict_orfs.targets[0], 
		blastp_sprot.targets[0],blastp_ur90_target,pfam.targets[0],signalp.targets[0],
		tmhmm.targets[0],dependencies)
	keg = tf.keg_task([annot])
	all_tasks = [busco,assembly_stats,gene_trans_map,blastx_sprot,rnammer,
				predict_orfs,signalp,blastp_sprot,tmhmm,pfam,annot,keg]
	if(cegma):
		all_tasks.append(cegma)
	if(blast_uniref90):
		all_tasks.append(blastp_ur90)
		all_tasks.append(blastx_ur90)
	return Supervisor(tasks=all_tasks,dependencies=dependency_set)


if(__name__=='__main__'):
	pass