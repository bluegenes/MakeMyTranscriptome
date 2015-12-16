import argparse
import os
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf


def cpumod(cpu, k): return int(round(float(cpu)/k))


def gen_annotation_supervisor(cpu, uniref90_flag, nr_flag, blast_flag, signalp_flag, tmhmm_flag, rnammer_flag, dependency_set):
    tasks = []
    annot_table_opts = {}
    def task_insert(task, name=None):
        tasks.append(task)
        if(name != None):
            annot_table_opts[name] = task.targets[0]

    gene_trans_map = tf.gene_trans_map_task([])
    task_insert(gene_trans_map, 'geneTransMap')
    predict_orfs = tf.predict_orfs_task(cpumod(cpu, 2), [])
    task_insert(predict_orfs, 'transdecoder')
    pfam = tf.pfam_task(cpumod(cpu, 2), [predict_orfs])
    task_insert(pfam, 'pfam')

    if(tmhmm_flag):
        tmhmm = tf.tmhmm_task([predict_orfs])
        task_insert(tmhmm, 'tmhmm')
    if(signalp_flag):
        signalp = tf.signalp_task([predict_orfs])
        task_insert(signalp, 'signalP')
    if(rnammer_flag):
        rnammer = tf.rnammer_task([])
        task_insert(rnammer, 'rnammer')
    if(blast_flag):
        blastx_sprot = tf.blastx_task(tf.PATH_SWISS_PROT, int(cpu/2), [])
        task_insert(blastx_sprot, 'spX')
        blastp_sprot = tf.blastp_task(tf.PATH_SWISS_PROT, int(cpu/2), [predict_orfs])
        task_insert(blastp_sprot, 'spP')
        if(uniref90_flag):
            blastx_ur90 = tf.blastx_task(tf.PATH_UNIREF90, cpumod(cpu, 2), [])
            task_insert(blastx_ur90, 'ur90X')
            blastp_ur90 = tf.blastp_task(tf.PATH_UNIREF90, cpumod(cpu, 2), [predict_orfs])
            task_insert(blastp_ur90, 'ur90P')
        if(nr_flag):
            blastx_nr = tf.blastx_task(tf.PATH_NR, cpumod(cpu, 2), [])
            task_insert(blastx_nr, 'nrX')
            blastp_nr = tf.blastp_task(tf.PATH_NR, cpumod(cpu, 2), [predict_orfs])
            task_insert(blastp_nr, 'nrP')
    else:
        dmnd_dependencies = []
        def dmnd_task_insert(task, name=None):
            dmnd_dependencies.append(task)
            task_insert(task, name)

        dmnd_xsprot = tf.diamondX_task(tf.PATH_SWISS_PROT, cpumod(cpu, 2), dmnd_dependencies[:])
        dmnd_task_insert(dmnd_xsprot)
        expand = tf.blast_augment_task(tf.PATH_SWISS_PROT, dmnd_xsprot.targets[0], [dmnd_xsprot])
        task_insert(expand, 'spX')
        dmnd_psprot = tf.diamondP_task(tf.PATH_SWISS_PROT, cpumod(cpu, 2), dmnd_dependencies+[predict_orfs])
        dmnd_task_insert(dmnd_psprot)
        expand = tf.blast_augment_task(tf.PATH_SWISS_PROT, dmnd_psprot.targets[0], [dmnd_psprot])
        task_insert(expand, 'spP')
        if(uniref90_flag):
            dmnd_xur90 = tf.diamondX_task(tf.PATH_UNIREF90, cpumod(cpu, 2), dmnd_dependencies[:])
            dmnd_task_insert(dmnd_xur90)
            expand = tf.blast_augment_task(tf.PATH_UNIREF90, dmnd_xur90.targets[0], [dmnd_xur90])
            task_insert(expand, 'ur90X')
            dmnd_pur90 = tf.diamondP_task(tf.PATH_UNIREF90, cpumod(cpu, 2), dmnd_dependencies+[predict_orfs])
            dmnd_task_insert(dmnd_pur90)
            expand = tf.blast_augment_task(tf.PATH_UNIREF90, dmnd_pur90.targets[0], [dmnd_pur90])
            task_insert(expand, 'ur90P')
        if(nr_flag):
            dmnd_xnr = tf.diamondX_task(tf.PATH_NR, cpumod(cpu, 2), dmnd_dependencies[:])
            dmnd_task_insert(dmnd_xnr)
            expand = tf.blast_augment_task(tf.PATH_NR, dmnd_xnr.targets[0], [dmnd_xnr])
            task_insert(expand, 'nrX')
            dmnd_pnr = tf.diamondP_task(tf.PATH_NR, cpumod(cpu, 2), dmnd_dependencies+[predict_orfs])
            dmnd_task_insert(dmnd_pnr)
            expand = tf.blast_augment_task(tf.PATH_NR, dmnd_pnr.targets[0], [dmnd_pnr])
            task_insert(expand, 'nrP')
    annot = tf.annot_table_task(annot_table_opts, tasks[:])
    tasks.append(annot)
    pipeplot = tf.pipeplot_task(annot.targets[0],[annot])
    tasks.append(pipeplot)
#    kegg = tf.kegg_task([annot])
#    tasks.append(kegg)
    return Supervisor(tasks=tasks,dependencies=dependency_set)


if(__name__=='__main__'):
    pass
