'''
'''

from tasks_v2 import Task
import os
from os.path import join, exists
import sys
import functions_general as fg
from external_tools import TOOLS_DICT
import re

''' static db variables '''
PATH_PFAM_DATABASE = '{0!s}/pfam/Pfam-A.hmm'.format(fg.PATH_DATABASES)
PATH_NR = join(fg.PATH_DATABASES, 'nr', 'nr')
PATH_SWISS_PROT = join(fg.PATH_DATABASES, 'uniprot_sprot', 'uniprot_sprot')
PATH_UNIREF90 = join(fg.PATH_DATABASES, 'uniref90', 'uniref90')
PATH_NOG_CATEGORIES = join(fg.PATH_DATABASES, 'nog_categories')

def gene_trans_map_task(path_assembly,out_dir,tasks):
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    trgs = ['{0!s}/{1!s}.gene_trans_map'.format(out_dir, assembly_name)]
    cmd = '{0!s} {1!s} > {2!s}'.format(fg.tool_path_check(TOOLS_DICT['trinity'].full_exe[1]),path_assembly,trgs[0])
    name = 'gene_trans_map_' + assembly_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def blast_task(blast_type, out_dir, path_query, path_db, cpu_cap, tasks):
    exe_index = 1
    if  blast_type == 'blastx':
        exe_index = 1
    elif blast_type == 'blastp':
        exe_index = 2
    assembly_name = os.path.basename(path_query).split('.')[0]
    db_name = os.path.basename(path_db).split('.')[0]
    trgs = ["{0!s}/{1!s}_{2!s}.{3!s}".format(out_dir, assembly_name, db_name, blast_type)] 
    cmd = ('{0!s} -query {1!s} -db {2!s} -num_threads {3!s} -max_target_seqs 1 '
            '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart '
            'send evalue bitscore stitle slen" -evalue 0.0001 > {4!s}'
            ).format(fg.tool_path_check(TOOLS_DICT['blast'].full_exe[exe_index]), path_query, path_db, cpu_cap, trgs[0])
    name = '{0!s}_{1!s}_{2!s}'.format(assembly_name, blast_type, db_name)
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,cpu=cpu_cap,stdout=out,stderr=err)

def diamond_task(blast_type, out_dir, path_query, ref, cpu_cap, tasks):
    ''' valid blast_types: "blastx", "blastp" ''' 
    base_ref = os.path.basename(ref)
    query_name = os.path.basename(path_query).split('.')[0]
    trgs = ['{0!s}/{1!s}_{2!s}.diamond_{3!s}'.format(out_dir, query_name, base_ref, blast_type)]
    pseudo_trgs = ['{0!s}/diamond_{1!s}_{2!s}'.format(out_dir, base_ref, blast_type)]
    cmd = ('{0!s} {1!s} --db {2!s} --query {3!s} --daa {4!s} --tmpdir {5!s} '
           '--max-target-seqs 20 --sensitive --threads {6!s} --evalue 0.001; {0!s} view '
           '--daa {4!s}.daa --out {7!s};').format(
           fg.tool_path_check(TOOLS_DICT['diamond'].full_exe[0]), blast_type, ref, path_query, pseudo_trgs[0], out_dir,
           cpu_cap, trgs[0])
    name = 'diamond_{0!s}_{1!s}_{2!s}'.format(blast_type, base_ref, query_name)
    out, err = fg.GEN_LOGS(name)
    return Task(command=cmd, dependencies=tasks, cpu=cpu_cap, targets=trgs, name=name, stdout=out, stderr=err)

def blast_augment_task(db, blast, tasks):
    id2name = db+'.stitle'
    trgs = ['{0!s}_ex'.format(blast)]
    cmd = 'python {0!s}/addStitleToBlastTab.py --db2Name {1!s} --blast {2!s} > {3!s}'.format(
           fg.PATH_SCRIPTS, id2name, blast, trgs[0])
    name = 'Blast_Augmentation_'+os.path.basename(blast)
    out, err = fg.GEN_LOGS(name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)

def rnammer_task(path_assembly, out_dir, tasks):
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    path_to_rnammer = os.path.dirname(TOOLS_DICT['rnammer'].folder_name)
    trgs = ['{0!s}/{1!s}.fasta.rnammer.gff'.format(out_dir,assembly_name)]
    cmd = ("cd {0!s}; {1!s} --transcriptome {2!s}  --path_to_rnammer {4!s} "
            "--org_type euk; cd -").format(out_dir,fg.tool_path_check(TOOLS_DICT['rnammer'].full_exe[0]),
            path_assembly,path_to_rnammer)
    name = 'rnammer_' + assembly_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def transdecoder_longorfs_task(path_assembly, path_transdecoder_output, cpu_cap, tasks):
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    longorf_outbase = os.path.join(path_transdecoder_output, assembly_name + '.fasta.transdecoder_dir') 
    trgs = ['{0!s}/longest_orfs.pep'.format(longorf_outbase),'{0!s}/longest_orfs.gff3'.format(longorf_outbase),
    '{0!s}/longest_orfs.cds'.format(longorf_outbase)]
    cmd = ("mkdir -p {0!s}; cd {0!s}; {1!s} -t {2!s}").format(path_transdecoder_output,
            fg.tool_path_check(TOOLS_DICT['transdecoder'].full_exe[0]),path_assembly,cpu_cap)
    name = 'TransDecoder_LongORFs_' + assembly_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap) 

def transdecoder_predict_orfs_task(path_assembly,path_transdecoder_output,tasks,pfam_input='',blastp_input=''):
    ''' Use transdecoder to predict ORF's from input fasta file. 
        Required for  downstream blastp, pfam, tmhmm, signalp.
        Targets: *transdecoder.*
    '''
    pfam,blastp,retain_blastp, retain_pfam = '','','',''
    if len(pfam_input) > 0:
        pfam = '--retain_pfam_hits ' + pfam_input
        retain_pfam = '_retain_pfam'
    if len(blastp_input) > 0:
        blastp = '--retain_blastp_hits ' + blastp_input
        retain_blastp = '_retain_blastp'
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    trgs = ['{0!s}/{1!s}.fasta.transdecoder.pep'.format(path_transdecoder_output,assembly_name)]
    cmd = ("mkdir -p {0!s}; cd {0!s}; {1!s} -t {2!s} {3!s} {4!s}").format(path_transdecoder_output,
            fg.tool_path_check(TOOLS_DICT['transdecoder'].full_exe[1]),path_assembly,pfam, blastp)
    name = 'TransDecoder_' + assembly_name + retain_pfam + retain_blastp
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def signalp_task(path_orfs,out_dir,tasks):
    out_name = os.path.basename(path_orfs).split('.')[0]
    trgs = ['{0!s}/{1!s}.signalp'.format(out_dir,out_name)]
    cmd = '{0!s} -f short -n {1!s} {2!s}'.format(fg.tool_path_check(TOOLS_DICT['signalp'].full_exe[0]),trgs[-1],path_orfs)
    name = 'signalp_' + out_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def tmhmm_task(path_orfs,out_dir,tasks):
    out_name = os.path.basename(path_orfs).split('.')[0]
    trgs = ['{0!s}/{1!s}.tmhmm'.format(out_dir, out_name)] 
    cmd = 'cd {0!s}; {1!s} --short < {2!s} > {3!s}'.format(out_dir,fg.tool_path_check(TOOLS_DICT['tmhmm'].full_exe[0]),path_orfs,trgs[0])
    name = 'tmhmm_' + out_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def pfam_task(path_orfs,out_dir, cpu_cap, tasks):
    out_name = os.path.basename(path_orfs).split('.')[0]
    trgs = ['{0!s}/{1!s}.pfam'.format(out_dir,out_name)]
    cmd = '{0!s} --cpu {1!s} --domtblout {2!s} {3!s} {4!s}'.format(
        fg.tool_path_check(TOOLS_DICT['hmmer'].full_exe[0]),cpu_cap,trgs[0],PATH_PFAM_DATABASE, path_orfs)
    name = 'pfam_' + out_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)

def annot_table_task(path_assembly,out_dir,opts, tasks):
    out_name = os.path.basename(path_assembly).split('.fa')[0]
    suffixes = ['annotation.txt','annotation_by_gene.txt']
    trgs = ['{0!s}/{1!s}_{2!s}'.format(os.path.dirname(path_assembly),out_name,sufx) for sufx in suffixes]
    cmd = (
        'python {0!s}/annot_table_main.py --fasta {1!s} --outfile {2!s}/{3!s} '
        '--ko2path {4!s}/orthology_pathway.list --sp2enzyme '
        '{4!s}/swiss_enzyme.list --enzyme2path {4!s}/enzyme_pathway.list '
        '--pfam2enzyme {4!s}/pfam_enzyme.list --go2path {4!s}/go_pathway.txt '
        '--nog2function {4!s}/allKOG_functional_info.txt '
        '--go2slim {4!s}/goslim_generic.obo --sp2ko {4!s}/idmapping.KO '
        '--sp2nog {4!s}/idmapping.eggNOG --sp2ortho {4!s}/idmapping.orthodb '
        '--sp2bioc {4!s}/idmapping.biocyc --sp2goentrez '
        '{4!s}/idmapping_selected.tab ').format(
        fg.PATH_SCRIPTS, path_assembly, os.path.dirname(path_assembly), out_name,
        fg.PATH_DATABASES)
    cmd += ' '.join(['--'+k+' '+opts[k] for k in opts])
    name = 'build_annotation_table_' + out_name
    out, err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def kegg_task(annotation_table, out_dir,tasks, kegg_map_id='ko01100'):
    ''' color pathways found in transcriptome on a given kegg map'''
    trgs = ['{0!s}/{1!s}.pdf'.format(out_dir,kegg_map_id),
            '{0!s}/{1!s}_KO.txt'.format(out_dir,kegg_map_id)]
    cmd = ('python {0!s}/color_pathways2.py --path {1!s} --transcriptomeKO {2!s}'
            ' --output {3!s}').format(fg.PATH_SCRIPTS,kegg_map_id,annotation_table,out_dir)
    name='draw_kegg_map_{0!s}_{1!s}'.format(os.path.basename(annotation_table), kegg_map_id)
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def pipeplot_task(annotation_table,out_dir,tasks):
    trgs = []
    cmd = 'mkdir -p {0!s}/plots ; cd {0!s}/plots ; python {1!s}/pipePlot.py -i {2!s} ;'.format(
            out_dir,fg.PATH_SCRIPTS,annotation_table)
    name = 'pipeplot'
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def assembly_to_bed_task(path_assembly,out_dir, tasks):
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    trgs = ['{0!s}/{1!s}.bed'.format(out_dir,assembly_name)]
    cmd = 'python {0!s}/fasta_to_bed_count_length.py {1!s} {2!s}'.format(
            fg.PATH_SCRIPTS,path_assembly,trgs[0])
    name = 'fasta_to_bed_' + assembly_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

