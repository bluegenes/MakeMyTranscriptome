import argparse
import os
from tasks_v2 import Supervisor, Task
import task_functions_v2 as tf


def cpumod(cpu, k):
    return int(round(float(cpu)/k))


def gen_annotation_supervisor(cpu, uniref90_flag, nr_flag, blast_flag, signalp_flag, tmhmm_flag, rnammer_flag, dependency_set):
    tasks = []
    gene_trans_map = tf.gene_trans_map_task([])
    tasks.append(gene_trans_map)
    predict_orfs = tf.predict_orfs_task(cpumod(cpu, 2), [])
    pep_path = predict_orfs.targets[0]
    tasks.append(predict_orfs)
    pfam = tf.pfam_task(pep_path, cpumod(cpu, 2), [predict_orfs])
    tasks.append(pfam)
    annot_table_opts = {'gene_trans_map': gene_trans_map.targets[0],
                        'predict_orfs': predict_orfs.targets[0],
                        'pfam': pfam.targets[0]}
    annot_table_dep = [gene_trans_map, predict_orfs, pfam]
    if(tmhmm_flag):
        tmhmm = tf.tmhmm_task(pep_path, [predict_orfs])
        tasks.append(tmhmm)
        annot_table_opts['tmhmm'] = tmhmm.targets[0]
        annot_table_dep.append(tmhmm)
    if(signalp_flag):
        signalp = tf.signalp_task(pep_path, [predict_orfs])
        tasks.append(signalp)
        annot_table_opts['signalp'] = signalp.targets[0]
        annot_table_dep.append(signalp)
    if(rnammer_flag):
        rnammer = tf.rnammer_task([])
        tasks.append(rnammer)
        annot_table_opts['rnammer'] = rnammer.targets[0]
        annot_table_dep.append(rnammer)
    if(blast_flag):
        blastx_sprot = tf.blastx_task(tf.PATH_SWISS_PROT, int(cpu/2), [])
        tasks.append(blastx_sprot)
        annot_table_opts['blastx_sp'] = blastx_sprot.targets[0]
        annot_table_dep.append(blastx_sprot)
        blastp_sprot = tf.blastp_task(pep_path, tf.PATH_SWISS_PROT, int(cpu/2), [predict_orfs])
        tasks.append(blastp_sprot)
        annot_table_opts['blastp_sp'] = blastp_sprot.targets[0]
        annot_table_dep.append(blastp_sprot)
        if(uniref90_flag):
            blastx_ur90 = tf.blastx_task(tf.PATH_UNIREF90, cpumod(cpu, 2), [])
            tasks.append(blastx_ur90)
            annot_table_opts['blastx_ur90'] = blastx_ur90.targets[0]
            annot_table_dep.append(blastx_ur90)
            blastp_ur90 = tf.blastp_task(pep_path, tf.PATH_UNIREF90, cpumod(cpu, 2), [predict_orfs])
            tasks.append(blastp_ur90)
            annot_table_opts['blastp_ur90'] = blastp_ur90.targets[0]
            annot_table_dep.append(blastp_ur90)
        if(nr_flag):
            blastx_nr = tf.blastx_task(tf.PATH_NR, cpumod(cpu, 2), [])
            tasks.append(blastx_nr)
            annot_table_dep['blastx_nr'] = blastx_nr.targets[0]
            annot_table_dep.append(blastx_nr)
            blastp_nr = tf.blastp_task(pep_path, tf.PATH_NR, cpumod(cpu, 2), [predict_orfs])
    else:
        dmnd_dependencies = []
        dmnd_xsprot = tf.diamond_task(tf.PATH_SWISS_PROT, 'blastx', cpumod(cpu, 2), [d for d in dmnd_dependencies])
        tasks.append(dmnd_xsprot)
        dmnd_dependencies.append(dmnd_xsprot)
        annot_table_opts['blastx_sp'] = dmnd_xsprot.targets[0]
        annot_table_dep.append(dmnd_xsprot)
        dmnd_psprot = tf.diamond_task(tf.PATH_SWISS_PROT, 'blastp', cpumod(cpu, 2), dmnd_dependencies+[predict_orfs], source=pep_path)
        dmnd_dependencies.append(dmnd_psprot)
        tasks.append(dmnd_psprot)
        annot_table_opts['blastx_sp'] = dmnd_psprot.targets[0]
        annot_table_dep.append(dmnd_psprot)
        if(uniref90_flag):
            dmnd_xur90 = tf.diamond_task(tf.PATH_UNIREF90, 'blastx', cpumod(cpu, 2), [d for d in dmnd_dependencies])
            dmnd_dependencies.append(dmnd_xur90)
            tasks.append(dmnd_xur90)
            annot_table_opts['blastx_ur90'] = dmnd_xur90.targets[0]
            annot_table_dep.append(dmnd_xur90)
            dmnd_pur90 = tf.diamond_task(tf.PATH_UNIREF90, 'blastp', cpumod(cpu, 2), dmnd_dependencies+[predict_orfs], source=pep_path)
            dmnd_dependencies.append(dmnd_pur90)
            tasks.append(dmnd_pur90)
            annot_table_opts['blastp_ur90'] = dmnd_pur90.targets[0]
        if(nr_flag):
            dmnd_xnr = tf.diamond_task(tf.PATH_NR, 'blastx', cpumod(cpu, 2), [d for d in dmnd_dependencies])
            tasks.append(dmnd_xnr)
            dmnd_dependencies.append(dmnd_xnr)
            annot_table_opts['blastx_nr'] = dmnd_xnr.targets[0]
            annot_table_dep.append(dmnd_xnr)
            dmnd_pnr = tf.diamond_task(tf.PATH_NR, 'blastp', cpumod(cpu, 2), dmnd_dependencies+[predict_orfs], source=pep_path)
            tasks.append(dmnd_pnr)
            annot_table_opts['blastp_nr'] = dmnd_pnr.targets[0]
            annot_table_dep.append(dmnd_pnr)
    annot = tf.annot_table_task(annot_table_opts, annot_table_dep)
    tasks.append(annot)
    pipeplot = tf.pipeplot_task(annot.targets[0],[annot])
    tasks.append(pipeplot)
    keg = tf.keg_task([annot])
    tasks.append(keg)
    return Supervisor(tasks=tasks,dependencies=dependency_set)


if(__name__=='__main__'):
    pass