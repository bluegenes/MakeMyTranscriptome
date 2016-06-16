from tasks_v2 import Task, Supervisor
import os
from functions_general import tool_path_check, gen_logs, make_dir_task, cp_task
from external_tools import TOOLS_DICT
import mmt_defaults as statics


# opc is the output_path_class object


def build_bowtie_task(opc, path_assembly, assembly_name, out_dir, tasks):
    trgs = ['{0!s}/{1!s}.1.bt2'.format(out_dir, assembly_name)]
    cmd = '{0!s} --offrate 1 -f {1!s} {2!s}/{3!s}'.format(
          tool_path_check(TOOLS_DICT['bowtie2'].full_exe[0]), path_assembly,
          out_dir, assembly_name)
    name = 'build_bowtie_' + assembly_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def bowtie2_unpaired_task(opc, bowtie2_index, out_dir, fastq, out_name, opt, cpu_cap, tasks):
    opts = ['-a -t --end-to-end', '-t --local']
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir, out_name)]
    cmd = ('{0!s} {1!s} -L {2!s} -N 1 --threads {3!s} -x {4!s} -U '
           '{5!s} | samtools view -Sb -').format(
           tool_path_check(TOOLS_DICT['bowtie2'].full_exe[1]),
           opts[opt], 22, cpu_cap, bowtie2_index, fastq)
    name = 'bowtie2_' + os.path.basename(bowtie2_index) + '_' + out_name
    out, err = gen_logs(opc.path_logs, name)
    out = trgs[0]
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap, shell=True)


def bowtie2_task(opc, bowtie2_index, out_dir, fastq1, fastq2, out_name, opt, cpu_cap, tasks):
    opts = ['-a -t --end-to-end', '-t --local']
    opts_name = ['express', 'intersectBed']
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir, out_name)]
    cmd = ('{0!s} {1!s} -L {2!s} -N 1 --maxins 800 --threads {3!s} -x {4!s} -1 '
           '{5!s} -2 {6!s} | samtools view -Sb -').format(
           tool_path_check(TOOLS_DICT['bowtie2'].full_exe[1]),
           opts[opt], 22, cpu_cap, bowtie2_index, fastq1, fastq2)
    name = 'bowtie2_' + os.path.basename(bowtie2_index) + '_' + out_name + '_' + opts_name[opt]
    out, err = gen_logs(opc.path_logs, name)
    out = trgs[0]
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap, shell=True)


def express_task(opc, bowtie2_index, assembly_path, out_dir, out_name, bam_input, tasks):
    trgs = ['{0!s}/{1!s}.xprs'.format(out_dir, out_name)]
    mkdir = make_dir_task(os.path.join(out_dir, out_name))
    pseudo_trgs = ['{0!s}/{1!s}/results.xprs'.format(out_dir, out_name)]
    cmd = ('{0!s} --output-dir {1!s}/{2!s} {3!s} {4!s}').format(
           tool_path_check(TOOLS_DICT['express'].full_exe[0]), out_dir,
           out_name, assembly_path, bam_input)
    name = 'express_' + os.path.basename(bowtie2_index) + '_' + out_name
    out, err = gen_logs(opc.path_logs, name)
    express = Task(command=cmd, dependencies=[make_dir_task], targets=pseudo_trgs, name=name, stdout=out, stderr=err)
    cp = cp_task(pseudo_trgs[0], trgs[0], [express])
    super_name = 'super_' + name
    return Supervisor(tasks=[mkdir, express, cp], dependencies=tasks, name=super_name)

def counts_to_table_task(opc, assembly_name, gene_trans_map, out_dir, count_files, out_name, flag, tasks):
    trgs = ['{0!s}/{1!s}.countsTable'.format(out_dir, out_name),
            '{0!s}/{1!s}_gene.countsTable'.format(out_dir, out_name)]
    count_str = ' '.join(['--counts {0!s}'.format(f) for f in count_files])
    cmd = ('python {0!s}/counts_to_table2.py --out {1!s} --inDir {2!s} '
           '--outDir {2!s} {3!s} {4!s} --geneTransMap {5!s}').format(
           statics.PATH_UTIL, out_name, out_dir, flag, count_str, gene_trans_map)
    if len(flag) > 1:
        name = '_' + flag.split('--')[1]
    else:
        name = '_intersectBed'
    name = 'counts_to_table_' + assembly_name + name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


# samtools is not currently included as a tools class instance
def sam_sort_task(opc, out_dir, bam_file, out_name, tasks):
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir, out_name)]
    cmd = 'samtools sort {0!s} {1!s}/{2!s}'.format(bam_file, out_dir, out_name)
    name = 'sam_sort_' + os.path.basename(bam_file) + '_' + out_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def intersect_bed_task(opc, out_dir, bam_file, bed_reference, output_name, tasks):
    trgs = ['{0!s}/{1!s}.bed'.format(out_dir, output_name)]
    # cmd = '{0!s} intersect -abam {1!s} -b {2!s} -wb -bed > {3!s}'.format(
    cmd = '{0!s} -abam {1!s} -b {2!s} -wb -bed'.format(
          tool_path_check(TOOLS_DICT['bedtools'].full_exe[0]),
          bam_file, bed_reference)
    name = 'intersect_bed_' + os.path.basename(bed_reference) + '_' + output_name
    out, err = gen_logs(opc.path_logs, name)
    out = trgs[0]
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def deseq2_task(opc, assembly_name, out_dir, counts_to_table_results, sample_info, basename, model, tasks):
    '''We don't check for R packages! need to figure this out..
    '''
    pseudo_model_temp = ''.join([c if(c != ' ') else '_' for c in model])
    trgs = ['{0!s}/deseq2_{1!s}_{2!s}/'.format(out_dir, basename, pseudo_model_temp)]
    cmd = 'Rscript {5!s}/deseq2.r --args {5!s} {0!s} {1!s} {2!s} {3!s} {4!s}'.format(
            counts_to_table_results, sample_info, out_dir, basename, model, statics.PATH_UTIL)
    # name = 'de_' + assembly_name + '_' + basename
    name = 'de_' + basename + '_' + os.path.basename(out_dir)
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def build_salmon_task(opc, path_assembly, assembly_name, out_dir, cpu_cap, tasks):
    trgs = ['{0!s}/{1!s}_salmon'.format(out_dir, assembly_name)]
    cmd = ('{0!s} index --transcripts {1!s} --index {2!s}/{3!s}_salmon '
           '--threads {4!s} --type quasi').format(
           tool_path_check(TOOLS_DICT['salmon'].full_exe[0]), path_assembly,
           out_dir, assembly_name, cpu_cap)
    name = 'build_salmon_' + assembly_name
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap)


def salmon_gene_map_task(opc, out_dir, assembly_name, gene_trans_map, tasks):
    ''' salmon requires gene_trans_map in reverse column order (transcript \\t gene \\n)'''
    trgs = ['{0!s}/{1!s}.trans_gene_map'.format(out_dir, assembly_name)]
    cmd = '''awk '{{ print $2 " " $1}}' {0!s}'''.format(gene_trans_map)
    name = 'salmon_gene_map_task_' + assembly_name
    out, err = gen_logs(opc.path_logs, name)
    out = trgs[0]
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def salmon_task(opc, index, left, right, out_name, gene_map, out_dir, cpu_cap, tasks):
    trgs = ['{0!s}/{1!s}_quant.sf'.format(out_dir, out_name),
            '{0!s}/{1!s}_quant.genes.sf'.format(out_dir, out_name)]
    pseudo_trgs = ['{0!s}/{1!s}/quant.sf'.format(out_dir, out_name),
                   '{0!s}/{1!s}/quant.genes.sf'.format(out_dir, out_name)]
    cmd = ('{0!s} quant -i {1!s} -l IU -1 {2!s} -2 {3!s} -o {4!s}/{5!s} '
           '--geneMap {6!s} -p {7!s} --extraSensitive').format(
    #cmd = '{0!s} quant -i {1!s} -l IU -1 {2!s} -2 {3!s} -o {4!s}/{5!s} 
    #--geneMap {6!s} -p {7!s} --extraSensitive --numBootstraps 30 --biasCorrect ; cp ' \
           tool_path_check(TOOLS_DICT['salmon'].full_exe[0]), index, left,
           right, out_dir, out_name, gene_map, cpu_cap)
    name = os.path.basename(index) + '_' + os.path.basename(left)
    out, err = gen_logs(opc.path_logs, name)
    salmon = Task(command=cmd, dependencies=[], targets=pseudo_trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap)
    cp1 = cp_task(pseudo_trgs[0], trgs[0], [salmon])
    cp2 = cp_task(pseudo_trgs[1], trgs[1], [salmon])
    super_name = 'super_' + name
    return Supervisor(tasks=[salmon, cp1, cp2], dependencies=tasks, name=super_name)



def salmon_unpaired_task(opc, index, unpaired, out_name, gene_map, out_dir, cpu_cap, tasks):
    trgs = ['{0!s}/{1!s}_quant.sf'.format(out_dir, out_name),
            '{0!s}/{1!s}_quant.genes.sf'.format(out_dir, out_name)]
    pseudo_trgs = ['{0!s}/{1!s}/quant.sf'.format(out_dir, out_name),
                   '{0!s}/{1!s}/quant.genes.sf'.format(out_dir, out_name)]
    cmd = ('{0!s} quant -i {1!s} -l U -r {2!s} -o {3!s}/{4!s} --geneMap {5!s} '
           '-p {6!s} --extraSensitive').format(
           tool_path_check(TOOLS_DICT['salmon'].full_exe[0]), index, unpaired,
           out_dir, out_name, gene_map, cpu_cap)
    name = 'salmon_unpaired_' + os.path.basename(index) + '_' + os.path.basename(unpaired)
    out, err = gen_logs(opc.path_logs, name)
    salmon = Task(command=cmd, dependencies=[], targets=pseudo_trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap)
    cp1 = cp_task(pseudo_trgs[0], trgs[0], [salmon])
    cp2 = cp_task(pseudo_trgs[1], trgs[1], [salmon])
    super_name = 'super_' + name
    return Supervisor(tasks=[salmon, cp1, cp2], dependencies=tasks, name=super_name)


def build_kallisto_task(opc, assembly_path, assembly_name, out_dir, tasks):
    # NO TARGETS
    trgs = []
    cmd = '{0!s} index -i {1!s}/{2!s}_kallisto {3!s}'.format(
          tool_path_check(TOOLS_DICT['kallisto'].full_exe[0]), out_dir,
          assembly_name, assembly_path)
    name = 'build_kallisto'
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def kallisto_task(opc, index,out_dir,out_name,left,right,tasks):
    # NO TARGETS
    trgs = []
    cmd = '{0!s} quant -i {1!s} -o {2!s}/{3!s} {4!s} {5!s}'.format(
          tool_path_check(TOOLS_DICT['kallisto'].full_exe[0]), index,
          out_dir, out_name, left, right)
    name = 'kallisto'
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)
