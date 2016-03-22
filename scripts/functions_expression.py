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
PATH_BUSCO_REFERENCE = join(fg.PATH_DATABASES, 'busco')
PATH_PFAM_DATABASE = '{0!s}/pfam/Pfam-A.hmm'.format(fg.PATH_DATABASES)
PATH_NR = join(fg.PATH_DATABASES, 'nr', 'nr')
PATH_SWISS_PROT = join(fg.PATH_DATABASES, 'uniprot_sprot', 'uniprot_sprot')
PATH_UNIREF90 = join(fg.PATH_DATABASES, 'uniref90', 'uniref90')
PATH_NOG_CATEGORIES = join(fg.PATH_DATABASES, 'nog_categories')

def build_bowtie_task(path_assembly, assembly_name, out_dir, tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.1.bt2'.format(out_dir,assembly_name)]  
    cmd = '{0!s} --offrate 1 -f {1!s} {2!s}/{3!s}'.format(
            fg.tool_path_check(TOOLS_DICT['bowtie2'].full_exe[0]), path_assembly, out_dir, assembly_name) 
    name = 'build_bowtie_' + assembly_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def bowtie2_unpaired_task(bowtie2_index,out_dir,fastq,out_name,opt,cpu_cap,tasks):
    '''
    '''
    opts = ['-a -t --end-to-end', '-t --local']
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir,out_name)]
    cmd = ('{0!s} {1!s} -L {2!s} -N 1 --threads {3!s} -x {4!s} -U '
            '{5!s} | samtools view -Sb - > {6!s} ').format(fg.tool_path_check(TOOLS_DICT['bowtie2'].full_exe[1]),
            opts[opt],22,cpu_cap,bowtie2_index,fastq,trgs[0])
    name = 'bowtie2_'+ os.path.basename(bowtie2_index) + '_' + out_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def bowtie2_task(bowtie2_index,out_dir,fastq1,fastq2,out_name,opt,cpu_cap,tasks):
    '''    
    '''
    opts = ['-a -t --end-to-end', '-t --local']
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir,out_name)]
    cmd = ('{0!s} {1!s} -L {2!s} -N 1 --maxins 800 --threads {3!s} -x {4!s} -1 '
            '{5!s} -2 {6!s} | samtools view -Sb - > {7!s} ').format(fg.tool_path_check(TOOLS_DICT['bowtie2'].full_exe[1]),
            opts[opt],22,cpu_cap,bowtie2_index,fastq1,fastq2,trgs[0])
    name = 'bowtie2_'+ os.path.basename(bowtie2_index) + '_' + out_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def express_task(assembly_path,out_dir,out_name,bam_input,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.xprs'.format(out_dir,out_name)]
    cmd = ('mkdir {1!s}/{2!s}; {0!s} --output-dir {1!s}/{2!s} {3!s} {4!s}; mv '
            '{1!s}/{2!s}/results.xprs {5!s}; rm -rf {1!s}/{2!s};').format(
            fg.tool_path_check(TOOLS_DICT['express'].full_exe[0]),out_dir,out_name,assembly_path,bam_input,trgs[0])
    name = 'express_'+ os.path.basename(bowtie2_index) + '_' + out_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def counts_to_table_task(assembly_name,gene_trans_map,out_dir,count_files,out_name,flag,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.countsTable'.format(out_dir,out_name),'{0!s}/{1!s}_byGene.countsTable'.format(out_dir,out_name)]
    count_str = ' '.join(['--counts {0!s}'.format(f) for f in count_files])
    cmd = ('python {0!s}/counts_to_table2.py --out {1!s} --inDir {2!s} '
            '--outDir {2!s} {3!s} {4!s} --geneTransMap {5!s}').format( 
            fg.PATH_SCRIPTS,out_name,out_dir,flag,count_str,gene_trans_map)
    name = 'counts_to_table_'+ assembly_name + '_' + flag.split('--')[1]
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

#samtools is not currently included as a tools class instance
def sam_sort_task(out_dir,bam_file,out_name,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir,out_name)]
    cmd = 'samtools sort {0!s} {1!s}/{2!s}'.format(bam_file,out_dir,out_name)
    name = 'sam_sort_'+ os.path.basename(bam_file) + '_' + out_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def intersect_bed_task(out_dir,bam_file,bed_reference,output_name,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.bed'.format(out_dir,output_name)]
    cmd = '{0!s} intersect -abam {1!s} -b {2!s} -wb -bed > {3!s}'.format(
        fg.tool_path_check(TOOLS_DICT['bedtools'].full_exe[0]),bam_file,bed_reference,trgs[0])
    name = 'intersect_bed_'+ bed_reference + '_' + output_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def deseq2_task(assembly_name,out_dir,counts_to_table_results,sample_info,basename,model,tasks):
    '''We don't check for R packages! need to figure this out..
    '''
    pseudo_model_temp = ''.join([c if c!=' ' else '_' for c in model])
    trgs = ['{0!s}/deseq2_{1!s}_{2!s}/'.format(out_dir,basename,pseudo_model_temp)]
    cmd = 'Rscript {5!s}/deseq2.r --args {5!s} {0!s} {1!s} {2!s} {3!s} {4!s}'.format(
            counts_to_table_results,sample_info,out_dir,basename,model,fg.PATH_SCRIPTS)
    name = 'de_' + assembly_name + '_' + basename 
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def build_salmon_task(path_assembly,assembly_name,out_dir,cpu_cap,tasks):
    trgs = ['{0!s}/{1!s}_salmon'.format(out_dir, assembly_name)] 
    cmd = '{0!s} index --transcripts {1!s} --index {2!s}/{3!s}_salmon --threads {4!s} --type quasi'.format(
        fg.tool_path_check(TOOLS_DICT['salmon'].full_exe[0]),path_assembly, out_dir, assembly_name, cpu_cap)
    name = 'build_salmon_' + assembly_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err, cpu=cpu_cap)

def salmon_gene_map_task(out_dir,assembly_name,gene_trans_map,tasks):
    ''' salmon requires gene_trans_map in reverse column order (transcript \\t gene \\n)'''
    trgs = ['{0!s}/{1!s}.trans_gene_map'.format(out_dir,assembly_name)] 
    cmd = 'join -t, -o 1.2,1.1 {0!s} {0!s} > {1!s}'.format(gene_trans_map, trgs[0]) 
    name = 'salmon_gene_map_task_' + assembly_name
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

def salmon_task(index,left,right,out_name,gene_map,out_dir,cpu_cap,tasks):
    trgs = ['{0!s}/{1!s}_quant.sf'.format(out_dir,out_name),'{0!s}/{1!s}_quant.genes.sf'.format(out_dir,out_name)]
    cmd = '{0!s} quant -i {1!s} -l IU -1 {2!s} -2 {3!s} -o {4!s}/{5!s} --geneMap {6!s} -p {7!s} --extraSensitive; cp ' \
        '{4!s}/{5!s}/quant.sf {4!s}/{5!s}_quant.sf; cp {4!s}/{5!s}/quant.genes.sf {4!s}/{5!s}_quant.genes.sf'.format(
	fg.tool_path_check(TOOLS_DICT['salmon'].full_exe[0]),index,left,right,out_dir,out_name,gene_map,cpu_cap)
    name = os.path.basename(index) + '_' + os.path.basename(left)
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def salmon_unpaired_task(index,unpaired,out_name,gene_map,out_dir,cpu_cap,tasks):
    trgs = ['{0!s}/{1!s}_quant.sf'.format(out_dir,out_name),'{0!s}/{1!s}_quant.genes.sf'.format(out_dir,out_name)]
    cmd = '{0!s} quant -i {1!s} -l U -r {2!s} -o {3!s}/{4!s} --geneMap {5!s} -p {6!s} --extraSensitive; cp ' \
        '{3!s}/{4!s}/quant.sf {3!s}/{4!s}_quant.sf; cp {3!s}/{4!s}/quant.genes.sf {3!s}/{4!s}_quant.genes.sf'.format(
        fg.tool_path_check(TOOLS_DICT['salmon'].full_exe[0]),index,unpaired,out_dir,out_name,gene_map,cpu_cap)
    name = 'salmon_unpaired_' + os.path.basename(index) + '_' + os.path.basename(unpaired)
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def build_kallisto_task(assembly_path, assembly_name,out_dir,tasks):
    trgs = []
    cmd = '{0!s} index -i {1!s}/{2!s}_kallisto {3!s}'.format(
            fg.tool_path_check(TOOLS_DICT['kallisto'].full_exe[0]),out_dir,assembly_name,assembly_path)
    name = 'build_kallisto'
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def kallisto_task(index,out_dir,out_name,left,right,tasks):
    trgs = []
    cmd = '{0!s} quant -i {1!s} -o {2!s}/{3!s} {4!s} {5!s}'.format(
            fg.tool_path_check(TOOLS_DICT['kallisto'].full_exe[0]),index,out_dir,out_name,left,right)
    name = 'kallisto'
    out,err = fg.GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)

