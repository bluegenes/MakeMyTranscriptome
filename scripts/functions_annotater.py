import os
from tasks_v2 import Task
import mmt_defaults as statics
from external_tools import TOOLS_DICT
from functions_general import gen_logs, tool_path_check


''' static db variables
PATH_PFAM_DATABASE = '{0!s}/pfam/Pfam-A.hmm'.format(fg.PATH_DATABASES)
PATH_NR = join(fg.PATH_DATABASES, 'nr', 'nr')
PATH_SWISS_PROT = join(fg.PATH_DATABASES, 'uniprot_sprot', 'uniprot_sprot')
PATH_UNIREF90 = join(fg.PATH_DATABASES, 'uniref90', 'uniref90')
PATH_NOG_CATEGORIES = join(fg.PATH_DATABASES, 'nog_categories')
'''

# opc is the output_path_class object

def gene_trans_map_task(opc, path_assembly, out_dir, tasks):
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    trgs = ['{0!s}/{1!s}.gene_trans_map'.format(out_dir, assembly_name)]
    cmd = 'perl {0!s} {1!s} > {2!s}'.format(
          tool_path_check(TOOLS_DICT['trinity'].full_exe[1]),
          path_assembly, trgs[0])
    name = 'gene_trans_map_' + assembly_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def blast_task(opc, blast_type, out_dir, path_query, path_db, cpu_cap, tasks):
    exe_index = 1
    if(blast_type == 'blastx'):
        exe_index = 1
    elif(blast_type == 'blastp'):
        exe_index = 2
    assembly_name = os.path.basename(path_query).split('.')[0]
    db_name = os.path.basename(path_db).split('.')[0]
    trgs = ["{0!s}/{1!s}_{2!s}.{3!s}".format(out_dir, assembly_name, db_name, blast_type)]
    cmd = ('{0!s} -query {1!s} -db {2!s} -num_threads {3!s} -max_target_seqs 1 '
           '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart '
           'send evalue bitscore stitle slen" -evalue 0.0001 > {4!s}'
           ).format(tool_path_check(TOOLS_DICT['blast'].full_exe[exe_index]),
                    path_query, path_db, cpu_cap, trgs[0])
    name = '{0!s}_{1!s}_{2!s}'.format(assembly_name, blast_type, db_name)
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, cpu=cpu_cap, stdout=out, stderr=err)


def diamond_task(opc, blast_type, out_dir, path_query, ref, cpu_cap, tasks):
    ''' valid blast_types: "blastx", "blastp" '''
    base_ref = os.path.basename(ref)
    query_name = os.path.basename(path_query).split('.')[0]
    trgs = ['{0!s}/{1!s}_{2!s}.diamond_{3!s}'.format(out_dir, query_name, base_ref, blast_type)]
    pseudo_trgs = ['{0!s}/diamond_{1!s}_{2!s}'.format(out_dir, base_ref, blast_type)]
    cmd = ('{0!s} {1!s} --db {2!s} --query {3!s} --daa {4!s} --tmpdir {5!s} '
           '--max-target-seqs 20 --sensitive --threads {6!s} --evalue 0.001; {0!s} view '
           '--daa {4!s}.daa --out {7!s};').format(
           tool_path_check(TOOLS_DICT['diamond'].full_exe[0]), blast_type, ref, path_query,
           pseudo_trgs[0], out_dir, cpu_cap, trgs[0])
    name = 'diamond_{0!s}_{1!s}_{2!s}'.format(blast_type, base_ref, query_name)
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, cpu=cpu_cap, targets=trgs, name=name, stdout=out, stderr=err)


def blast_augment_task(opc, db, blast, tasks):
    id2name = db+'.stitle'
    trgs = ['{0!s}_ex'.format(blast)]
    cmd = 'python {0!s}/addStitleToBlastTab.py --db2Name {1!s} --blast {2!s} > {3!s}'.format(
           statics.PATH_UTIL, id2name, blast, trgs[0])
    name = 'Blast_Augmentation_'+os.path.basename(blast)
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def rnammer_task(opc, path_assembly, out_dir, tasks):
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    path_to_rnammer = os.path.dirname(TOOLS_DICT['rnammer'].folder_name)
    trgs = ['{0!s}/{1!s}.fasta.rnammer.gff'.format(out_dir, assembly_name)]
    cmd = ("cd {0!s}; {1!s} --transcriptome {2!s}  --path_to_rnammer {4!s} "
           "--org_type euk; cd -").format(
           out_dir, tool_path_check(TOOLS_DICT['rnammer'].full_exe[0]),
           path_assembly, path_to_rnammer)
    name = 'rnammer_' + assembly_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def transdecoder_longorfs_task(opc, path_assembly, path_transdecoder_output, cpu_cap, tasks):
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    longorf_outbase = os.path.join(path_transdecoder_output, opc.assembly_name + '.fasta.transdecoder_dir')
    trgs = ['{0!s}/longest_orfs.pep'.format(longorf_outbase),
            '{0!s}/longest_orfs.gff3'.format(longorf_outbase),
            '{0!s}/longest_orfs.cds'.format(longorf_outbase)]
    cmd = ("mkdir -p {0!s}; cd {0!s}; {1!s} -t {2!s}").format(path_transdecoder_output,
            tool_path_check(TOOLS_DICT['transdecoder'].full_exe[0]), path_assembly, cpu_cap)
    name = 'TransDecoder_LongORFs_' + assembly_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap) 


def transdecoder_predict_orfs_task(opc, path_assembly, path_transdecoder_output, tasks, pfam_input='', blastp_input=''):
    ''' Use transdecoder to predict ORF's from input fasta file. 
        Required for  downstream blastp, pfam, tmhmm, signalp.
        Targets: *transdecoder.*
    '''
    pfam, blastp, retain_blastp, retain_pfam = '', '', '', ''
    if len(pfam_input) > 0:
        pfam = '--retain_pfam_hits ' + pfam_input
        retain_pfam = '_retain_pfam'
    if len(blastp_input) > 0:
        blastp = '--retain_blastp_hits ' + blastp_input
        retain_blastp = '_retain_blastp'
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    base_targ = '{0!s}/{1!s}.fasta.transdecoder'.format(path_transdecoder_output, opc.assembly_name)
    trgs = [base_targ+'.pep', base_targ+'.bed', base_targ+'.gff3']
    cmd = "mkdir -p {0!s}; cd {0!s}; {1!s} -t {2!s} {3!s} {4!s}".format(
          path_transdecoder_output, tool_path_check(TOOLS_DICT['transdecoder'].full_exe[1]),
          path_assembly, pfam, blastp)
    name = 'TransDecoder_Predict_' + assembly_name + retain_pfam + retain_blastp
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def signalp_task(opc, path_orfs, out_dir, tasks):
    out_name = os.path.basename(path_orfs).split('.')[0]
    trgs = ['{0!s}/{1!s}.signalp'.format(out_dir, out_name)]
    cmd = '{0!s} -f short -n {1!s} {2!s}'.format(tool_path_check(TOOLS_DICT['signalp'].full_exe[0]), trgs[-1], path_orfs)
    name = 'signalp_' + out_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def tmhmm_task(opc, path_orfs, out_dir, tasks):
    out_name = os.path.basename(path_orfs).split('.')[0]
    trgs = ['{0!s}/{1!s}.tmhmm'.format(out_dir, out_name)]
    cmd = 'cd {0!s}; {1!s} --short < {2!s} > {3!s}'.format(
          out_dir, tool_path_check(TOOLS_DICT['tmhmm'].full_exe[0]),
          path_orfs, trgs[0])
    name = 'tmhmm_' + out_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def pfam_task(opc, dbs, path_orfs, out_dir, cpu_cap, tasks):
    out_name = os.path.basename(path_orfs).split('.')[0]
    trgs = ['{0!s}/{1!s}.pfam'.format(out_dir, out_name)]
    cmd = '{0!s} --cpu {1!s} --domtblout {2!s} {3!s} {4!s}'.format(
          tool_path_check(TOOLS_DICT['hmmer'].full_exe[2]), cpu_cap,
          trgs[0], dbs['pfam'].call_path, path_orfs)
    name = 'pfam_' + out_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap)


def pfam_seq_task(opc, dbs, path_orfs, out_dir, cpu_cap, tasks):
    out_name = os.path.basename(path_orfs).split('.')[0]
    trgs = ['{0!s}/{1!s}.pfam_tblout'.format(out_dir, out_name)]
    cmd = '{0!s} --cpu {1!s} --tblout {2!s} {3!s} {4!s}'.format(
          tool_path_check(TOOLS_DICT['hmmer'].full_exe[0]), cpu_cap,
          trgs[0], dbs['pfam'].call_path, path_orfs)
    name = 'pfam_tblout_' + out_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap)


def annot_table_task(opc, dbs, path_assembly, out_dir, opts, tasks):
    out_name = os.path.basename(path_assembly).split('.fa')[0]
    suffixes = ['annotation.txt']  # ,'annotation_by_gene.txt']
    base_out_name = '{0!s}/{1!s}'.format(os.path.dirname(path_assembly), out_name)
    trgs = ['{0!s}_{1!s}'.format(base_out_name, sufx) for sufx in suffixes]
    grab_db = lambda s: dbs[s].call_path
    cmd = (
        'python {0!s}/annot_table_pandas.py --fasta {1!s} --outfile {2!s} '
        '--ko2path {14!s} --sp2enzyme {3!s} '
        '--enzyme2path {4!s} --pfam2enzyme {5!s} --go2path {6!s} '
        '--nog2function {7!s} --go2slim {8!s} --sp2ko {9!s} --sp2nog {10!s}'
        ' --sp2ortho {11!s} --sp2bioc {12!s} --sp2goentrez {13!s} ').format(
        statics.PATH_UTIL, path_assembly, base_out_name, grab_db('swiss_enzyme'),
        grab_db('enzyme_pathway'), grab_db('pfam_enzyme'), grab_db('go_pathway'),
        grab_db('nog_categories'), grab_db('goslim_generic'), grab_db('id_mapping_ko'),
        grab_db('id_mapping_eggnog'), grab_db('id_mapping_orthodb'),
        grab_db('id_mapping_biocyc'), grab_db('id_mapping_selected'), grab_db('orthology_pathway'))
    cmd += ' '.join(['--'+k+' '+opts[k] for k in opts])
    name = 'build_annotation_table_' + out_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def gff3_task(opc, path_assembly, out_path, opts, tasks):
    trgs = [out_path]
    cmd = ('python {0!s}/annot_table_gff3.py --fasta {1!s} --outfile {2!s} '
           ).format(statics.PATH_UTIL, path_assembly, out_path)
    cmd += ' '.join(['--'+k+' '+opts[k] for k in opts])
    name = 'build_gff3_' + os.path.basename(out_path)
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def kegg_task(opc, annotation_table, out_dir, tasks, kegg_map_id='ko01100'):
    kegg_dir = '{0!s}/kegg_maps'.format(out_dir)
    trgs = ['{0!s}/{1!s}.pdf'.format(kegg_dir, kegg_map_id),
            '{0!s}/{1!s}_KO.txt'.format(kegg_dir, kegg_map_id)]
    cmd = ('mkdir -p {3!s} ; cd {3!s} ; python {0!s}/color_pathways2.py --path {1!s} '
           ' --transcriptomeKO {2!s} --output {3!s}').format(
           statics.PATH_UTIL, kegg_map_id, annotation_table, kegg_dir)
    name = 'draw_kegg_map_{0!s}_{1!s}'.format(os.path.basename(annotation_table), kegg_map_id)
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def pipeplot_task(opc, dbs, annotation_table, out_dir, tasks):
    trgs = ['{0!s}/plots/cogMultiple.png'.format(out_dir)]
    # pipeplot no targets
    cmd = 'mkdir -p {0!s}/plots ; cd {0!s}/plots ; python {1!s}/pipePlot.py -i {2!s} --nog_categories {3!s};'.format(
            out_dir, statics.PATH_UTIL, annotation_table, dbs['nog_categories'].call_path)
    name = 'pipeplot'
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def assembly_to_bed_task(opc, path_assembly, out_dir, tasks):
    assembly_name = os.path.basename(path_assembly).split('.fa')[0]
    trgs = ['{0!s}/{1!s}.bed'.format(out_dir, assembly_name)]
    cmd = 'python {0!s}/fasta_to_bed_count_length.py {1!s} {2!s}'.format(
            statics.PATH_UTIL, path_assembly, trgs[0])
    name = 'fasta_to_bed_' + assembly_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)
