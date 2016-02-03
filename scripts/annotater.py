import argparse
import os
from tasks_v2 import Supervisor, Task
#import task_functions_v2 as tf
import functions_general as fg
import functions_annotater as fan
import functions_databases as fd


def cpumod(cpu, k): return int(round(float(cpu)/k))


def gen_annotation_supervisor(cpu, uniref90_flag, nr_flag, blast_flag, signalp_flag, tmhmm_flag, rnammer_flag, dependency_set, path_assembly=fg.GEN_PATH_ASSEMBLY(), assembly_name=fg.NAME_ASSEMBLY,out_dir=fg.GEN_PATH_ANNOTATION_FILES()):
    tasks = []
    annot_table_opts = {}
    def task_insert(task, name=None):
        tasks.append(task)
        if(name != None):
            annot_table_opts[name] = task.targets[0]
    gene_trans_map = fan.gene_trans_map_task(path_assembly,out_dir,[])
    task_insert(gene_trans_map, 'geneTransMap')
    predict_orfs = fan.predict_orfs_task(path_assembly, out_dir, cpumod(cpu, 2), [])
    task_insert(predict_orfs, 'transdecoder')
    pfam = fan.pfam_task(predict_orfs.targets[0], out_dir,cpumod(cpu, 2), [predict_orfs])
#    task_insert(pfam, 'pfam') # having trouble with pfam parsing errors
    tasks.append(pfam)
    if(blast_flag):
        blastx_sprot = fan.blast_task('blastx', out_dir, path_assembly, fd.PATH_SWISS_PROT, int(cpu/2), [])
        task_insert(blastx_sprot, 'spX')
        blastp_sprot = fan.blast_task('blastp',out_dir, predict_orfs.targets[0],fd.PATH_SWISS_PROT, int(cpu/2), [predict_orfs])
        task_insert(blastp_sprot, 'spP')

        if(uniref90_flag):
            blastx_ur90 = fan.blast_task('blastx',out_dir, path_assembly,fd.PATH_UNIREF90, cpumod(cpu, 2), [])
            task_insert(blastx_ur90, 'ur90X')
            blastp_ur90 = fan.blast_task('blastp',out_dir, predict_orfs.targets[0],fd.PATH_UNIREF90, cpumod(cpu, 2), [predict_orfs])
            task_insert(blastp_ur90, 'ur90P')
        if(nr_flag):
            blastx_nr = fan.blast_task('blastx',out_dir, path_assembly,fd.PATH_NR, cpumod(cpu, 2), [])
            task_insert(blastx_nr, 'nrX')
            blastp_nr = fan.blast_task('blastp',out_dir, predict_orfs.targets[0],fd.PATH_NR, cpumod(cpu, 2), [predict_orfs])
            task_insert(blastp_nr, 'nrP')
    else:
        dmnd_dependencies = []
        def dmnd_task_insert(task, name=None):
            dmnd_dependencies.append(task)
            task_insert(task, name)
        dmnd_xsprot = fan.diamond_task('blastx',out_dir, path_assembly,fd.PATH_SWISS_PROT, cpumod(cpu, 2), dmnd_dependencies[:])
        dmnd_task_insert(dmnd_xsprot)
        expand = fan.blast_augment_task(fd.PATH_SWISS_PROT, dmnd_xsprot.targets[0], [dmnd_xsprot])
        task_insert(expand, 'spX')
        dmnd_psprot = fan.diamond_task('blastp',out_dir, predict_orfs.targets[0],fd.PATH_SWISS_PROT, cpumod(cpu, 2), dmnd_dependencies+[predict_orfs])
        dmnd_task_insert(dmnd_psprot)
        expand = fan.blast_augment_task(fd.PATH_SWISS_PROT, dmnd_psprot.targets[0], [dmnd_psprot])
        task_insert(expand, 'spP')
        if(uniref90_flag):
            dmnd_xur90 = fan.diamond_task('blastx',out_dir,path_assembly,fd.PATH_UNIREF90, cpumod(cpu, 2), dmnd_dependencies[:])
            dmnd_task_insert(dmnd_xur90)
            expand = fan.blast_augment_task(fd.PATH_UNIREF90, dmnd_xur90.targets[0], [dmnd_xur90])
            task_insert(expand, 'ur90X')
            dmnd_pur90 = fan.diamond_task('blastp',out_dir, predict_orfs.targets[0],fd.PATH_UNIREF90, cpumod(cpu, 2), dmnd_dependencies+[predict_orfs])
            dmnd_task_insert(dmnd_pur90)
            expand = fan.blast_augment_task(fd.PATH_UNIREF90, dmnd_pur90.targets[0], [dmnd_pur90])
            task_insert(expand, 'ur90P')
        if(nr_flag):
            dmnd_xnr = fan.diamond_task('blastx',out_dir, path_assembly,fd.PATH_NR, cpumod(cpu, 2), dmnd_dependencies[:])
            dmnd_task_insert(dmnd_xnr)
            expand = fan.blast_augment_task(fd.PATH_NR, dmnd_xnr.targets[0], [dmnd_xnr])
            task_insert(expand, 'nrX')
            dmnd_pnr = fan.diamondP_task('blastp',out_dir, predict_orfs.targets[0], fd.PATH_NR, cpumod(cpu, 2), dmnd_dependencies+[predict_orfs])
            dmnd_task_insert(dmnd_pnr)
            expand = fan.blast_augment_task(fd.PATH_NR, dmnd_pnr.targets[0], [dmnd_pnr])
            task_insert(expand, 'nrP')
    if(tmhmm_flag):
        tmhmm = fan.tmhmm_task(predict_orfs.targets[0], out_dir,[predict_orfs])
        task_insert(tmhmm, 'tmhmm')
    if(signalp_flag):
        signalp = fan.signalp_task(predict_orfs.targets[0], out_dir, [predict_orfs])
        task_insert(signalp, 'signalP')
    if(rnammer_flag):
        rnammer = fan.rnammer_task(path_assembly,[])
        task_insert(rnammer, 'rnammer')
    
    annot = fan.annot_table_task(path_assembly,out_dir,annot_table_opts, tasks[:])
    tasks.append(annot)
    pipeplot = fan.pipeplot_task(annot.targets[0],out_dir,[annot])
    tasks.append(pipeplot)
#    kegg = fan.kegg_task([annot])
#    tasks.append(kegg)
    return Supervisor(tasks=tasks,dependencies=dependency_set)


if(__name__=='__main__'):
    pass
