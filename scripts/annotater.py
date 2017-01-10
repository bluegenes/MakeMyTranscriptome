import argparse
import os
from tasks_v2 import Supervisor, Task
#import task_functions_v2 as tf
import functions_general as fg
import functions_annotater as fan
import functions_expression as fex
# import functions_databases as fd


def cpumod(cpu, k): return int(round(float(cpu)/k))


def gen_annotation_supervisor(opc, dbs, cpu, uniref90_flag, nr_flag, blast_flag, signalp_flag, tmhmm_flag, rnammer_flag, dependency_set, gene_trans_map, path_assembly, assembly_name,out_dir, improve_orfs=False):
    tasks = []
    annot_table_opts = {'geneTransMap':gene_trans_map}
    gff3_dependencies = []
    gff3_opts = {}
    def task_insert(task, name=None, index=0, gff3_flag=False):
        tasks.append(task)
        if(name != None):
            annot_table_opts[name] = task.targets[index]
        if(gff3_flag):
            gff3_dependencies.append(task)
            gff3_opts[name] = task.targets[index]
    annot_table_opts['geneTransMap'] = gene_trans_map
    transd_dir = os.path.join(out_dir,'transdecoder')
    longorfs = fan.transdecoder_longorfs_task(opc, path_assembly,  transd_dir, cpumod(cpu, 2), [])
    tasks.append(longorfs)
    if improve_orfs:
        blastp_transd = fan.blast_task(opc, 'blastp',  transd_dir, longorfs.targets[0],dbs['uniprot_sprot'].call_path, int(cpu/2), [longorfs])
        pfam_transd = fan.pfam_task(opc,dbs,longorfs.targets[0], transd_dir,cpumod(cpu,2), [longorfs])
        tasks.extend([blastp_transd,pfam_transd]) 
        predict_orfs=fan.transdecoder_predict_orfs_task(opc, path_assembly,transd_dir,[longorfs,pfam_transd,blastp_transd],pfam_transd.targets[0],blastp_transd.targets[0])
    else:
        predict_orfs = fan.transdecoder_predict_orfs_task(opc, path_assembly,transd_dir,[longorfs])
    gff3_dependencies.append(predict_orfs)
    gff3_opts['transdecoder_gff3'] = predict_orfs.targets[2]
    task_insert(predict_orfs, 'transdecoder', 1)
    pfam = fan.pfam_task(opc, dbs, predict_orfs.targets[0], out_dir,cpumod(cpu, 4), [predict_orfs])
    #pfam = fan.pfam_task(predict_orfs.targets[0], out_dir,cpu, [predict_orfs])
    task_insert(pfam, 'pfam', gff3_flag=True) 
    if(blast_flag):
        blastx_sprot = fan.blast_task(opc, 'blastx', out_dir, path_assembly, dbs['uniprot_sprot'].call_path, cpumod(cpu, 2), [])
        task_insert(blastx_sprot, 'spX', gff3_flag=True)
        blastp_sprot = fan.blast_task(opc, 'blastp',out_dir, predict_orfs.targets[0], dbs['uniprot_sprot'].call_path, cpumod(cpu, 2), [predict_orfs])
        task_insert(blastp_sprot, 'spP', gff3_flag=True)

        if(uniref90_flag):
            blastx_ur90 = fan.blast_task(opc, 'blastx',out_dir, path_assembly,dbs['uniref90'].call_path, cpumod(cpu, 2), [])
            task_insert(blastx_ur90, 'ur90X', gff3_flag=True)
            blastp_ur90 = fan.blast_task(opc, 'blastp',out_dir, predict_orfs.targets[0],dbs['uniref90'].call_path, cpumod(cpu, 2), [predict_orfs])
            task_insert(blastp_ur90, 'ur90P', gff3_flag=True)
        if(nr_flag):
            blastx_nr = fan.blast_task(opc, 'blastx',out_dir, path_assembly,dbs['nr'].call_path, cpumod(cpu, 2), [])
            task_insert(blastx_nr, 'nrX', gff3_flag=True)
            blastp_nr = fan.blast_task(opc, 'blastp',out_dir, predict_orfs.targets[0],dbs['nr'].call_path, cpumod(cpu, 2), [predict_orfs])
            task_insert(blastp_nr, 'nrP', gff3_flag=True)
    else:
        dmnd_dependencies = []
        def dmnd_task_insert(task, name=None):
            dmnd_dependencies.append(task)
            task_insert(task, name)
        dmnd_xsprot = fan.diamond_task(opc, 'blastx',out_dir, path_assembly,dbs['uniprot_sprot'].call_path, cpumod(cpu, 2), dmnd_dependencies[:])
        dmnd_task_insert(dmnd_xsprot)
        expand = fan.blast_augment_task(opc, dbs['uniprot_sprot'].call_path, dmnd_xsprot.targets[0], [dmnd_xsprot])
        task_insert(expand, 'spX', gff3_flag=True)
        dmnd_psprot = fan.diamond_task(opc, 'blastp',out_dir, predict_orfs.targets[0],dbs['uniprot_sprot'].call_path, cpumod(cpu, 2), dmnd_dependencies+[predict_orfs])
        dmnd_task_insert(dmnd_psprot)
        expand = fan.blast_augment_task(opc, dbs['uniprot_sprot'].call_path, dmnd_psprot.targets[0], [dmnd_psprot])
        task_insert(expand, 'spP', gff3_flag=True)
        if(uniref90_flag):
            dmnd_xur90 = fan.diamond_task(opc, 'blastx',out_dir,path_assembly,dbs['uniref90'].call_path, cpumod(cpu, 2), dmnd_dependencies[:])
            dmnd_task_insert(dmnd_xur90)
            expand = fan.blast_augment_task(opc, dbs['uniref90'].call_path, dmnd_xur90.targets[0], [dmnd_xur90])
            task_insert(expand, 'ur90X', gff3_flag=True)
            dmnd_pur90 = fan.diamond_task(opc, 'blastp',out_dir, predict_orfs.targets[0],dbs['uniref90'].call_path, cpumod(cpu, 2), dmnd_dependencies+[predict_orfs])
            dmnd_task_insert(dmnd_pur90)
            expand = fan.blast_augment_task(opc, dbs['uniref90'].call_path, dmnd_pur90.targets[0], [dmnd_pur90])
            task_insert(expand, 'ur90P', gff3_flag=True)
        if(nr_flag):
            dmnd_xnr = fan.diamond_task(opc, 'blastx',out_dir, path_assembly,dbs['nr'].call_path, cpumod(cpu, 2), dmnd_dependencies[:])
            dmnd_task_insert(dmnd_xnr)
            expand = fan.blast_augment_task(opc, dbs['nr'].call_path, dmnd_xnr.targets[0], [dmnd_xnr])
            task_insert(expand, 'nrX', gff3_flag=True)
            dmnd_pnr = fan.diamond_task(opc, 'blastp',out_dir, predict_orfs.targets[0], dbs['nr'].call_path, cpumod(cpu, 2), dmnd_dependencies+[predict_orfs])
            dmnd_task_insert(dmnd_pnr)
            expand = fan.blast_augment_task(opc, dbs['nr'].call_path, dmnd_pnr.targets[0], [dmnd_pnr])
            task_insert(expand, 'nrP', gff3_flag=True)
    if(tmhmm_flag):
        tmhmm = fan.tmhmm_task(opc,predict_orfs.targets[0], out_dir,[predict_orfs])
        task_insert(tmhmm, 'tmhmm')
    if(signalp_flag):
        signalp = fan.signalp_task(opc,predict_orfs.targets[0], out_dir, [predict_orfs])
        task_insert(signalp, 'signalP')
    # need more intelligent annot table -- if pfam fails, for example, we can still generate an annot table
    annot = fan.annot_table_task(opc, dbs, path_assembly,out_dir,annot_table_opts, tasks[:])
    tasks.append(annot)
    gff3_output = os.path.join(opc.path_dir, opc.assembly_name+'.gff3')
    gff3 = fan.gff3_task(opc, path_assembly, gff3_output, gff3_opts, gff3_dependencies)
    tasks.append(gff3)
    pipeplot = fan.pipeplot_task(opc, dbs, annot.targets[0],out_dir,[annot])
    tasks.append(pipeplot)
    kegg = fan.kegg_task(opc, annot.targets[0],out_dir, [annot])
    tasks.append(kegg)
    return Supervisor(tasks=tasks,dependencies=dependency_set)


if(__name__=='__main__'):
    pass
