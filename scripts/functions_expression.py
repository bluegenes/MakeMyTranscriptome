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

def cegma_task(out_dir,assembly,cpu_cap, tasks):
    '''    Defines the cegma task. Uses PATH_DIR, PATH_CEGMA, NAME_ASSEMBLY.
        Params :
            cpu_cap - number of threads to be used by cegma
            tasks - a list of tasks that this task is dependant on (trinity_task)
    '''
    assembly_name = os.path.basename(assembly).split('.fa')[0]
    trgs = ['{0!s}/{1!s}.completeness_report'.format(out_dir,assembly_name)]
    cmd = '{0!s} -g {1!s} -v -o {3!s}/{2!s} -T {4!s}'.format(PATH_CEGMA,
            assembly,assembly_name,out_dir,cpu_cap)
    name = 'cegma'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,cpu=cpu_cap,stdout=out,stderr=err)


def busco_task(assembly_path, assembly_name, out_dir,reference_name, cpu_cap, tasks):
    ''' Defines the busco task. Uses PATH_DIR, PATH_BUSCO, PATH_BUSCO_REFERENCE
        Params :
            reference_name - Name of the reference file to be used by busco
            cpu_cap - the cpu limit to be gicen to busco.
            tasks - a list of tasks that this task is dependant on.
    '''
    trgs = ['{0!s}/run_busco_{1!s}_{2!s}'.format(out_dir,assembly_name,reference_name)]
    cmd = ('cd {0!s}; /matta1/biotools/anaconda/envs/py3k/bin/python {1!s} '
            '-o busco_{3!s}_{2!s} -in {4!s} -l {5!s}/{2!s}_buscos/{2!s} -m trans -f -c {6!s}'
            ).format(out_dir,tool_path_check(TOOLS_DICT['busco_plant'].full_exe[0]),reference_name,assembly_name,assembly_path,
            PATH_BUSCO_REFERENCE,cpu_cap)
    name = 'busco_'+ reference_name + '_' + assembly_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,cpu=cpu_cap,stdout=out,stderr=err)


def transrate_dep_generator(transrate_task, lefts, rights, singles, reference, assembly_path, cpu_cap, transrate_dir, other_dependencies):

    def ret():
        for t in other_dependencies:
            try:
                if( not d.finished()):
                    return False
            except Task.ExitCodeException:
                return False
        assembly_files = sorted(os.listdir(GEN_PATH_ASSEMBLY_FILES()))
        assembly_files = [os.path.join(GEN_PATH_ASSEMBLY_FILES(), f) for f in assembly_files]
        new_lefts = [[g for g in assembly_files if(os.path.basename(f) in g)] for f in lefts]
        new_lefts = [k[0] for k in new_lefts if(len(k) > 0)]
        new_rights = [[g for g in assembly_files if(os.path.basename(f) in g)] for f in rights]
        new_rights = [k[0] for k in new_rights if(len(k) > 0)]
        new_singles = [[g for g in assembly_files if(os.path.basename(f) in g)] for f in singles]
        new_singles = [k[0] for k in new_singles if(len(k) > 0)]
        if(len(new_lefts) == len(lefts) and len(new_rights) == len(rights) and len(new_singles) == len(singles)
            and len(new_lefts)+len(new_singles) != 0):
            new_lefts = ','.join(new_lefts+new_singles)
            new_rights = ','.join(new_rights) 
            new_lefts = '--left '+new_lefts if(len(new_lefts) > 0) else ''
            new_rights = '--right '+new_rights if(len(new_rights) > 0) else ''
            cmd = '{0!s} --assembly {1!s} {2!s} {3!s} --threads {4!s} {5!s} --output {6!s}'.format(
                   tool_path_check(TOOLS_DICT['transrate'].full_exe[0]), assembly_path, new_lefts,
                   new_rights, cpu_cap, reference, transrate_dir)
            transrate_task.command = cmd
        else:
            print('Unable to match input files with trimmed output. Continuing transrate using input files instead.')
        return True
    return ret


def transrate_task(assembly_path, assembly_name,lefts, rights, singles, out_dir, transrate_dir, cpu_cap, tasks, reference = ''): #, cpu_cap, tasks):
    trgs = ['{0!s}/assemblies.csv'.format(transrate_dir),'{0!s}/{1!s}/good.{1!s}.fasta'.format(transrate_dir,assembly_name),'{0!s}/{1!s}/{1!s}.fasta_quant.sf'.format(transrate_dir,assembly_name)]
    orig_lefts = lefts
    orig_rights = rights
    orig_singles = singles
    lefts = ','.join(lefts+singles)
    rights = ','.join(rights) 
    lefts = '--left '+lefts if(len(lefts) > 0) else ''
    rights = '--right '+rights if(len(rights) > 0) else ''
    reference = '--reference ' + reference if(reference != '') else ''
    cmd = '{0!s} --assembly {1!s} {2!s} {3!s} --threads {4!s} {5!s} --output {6!s}'.format(
           tool_path_check(TOOLS_DICT['transrate'].full_exe[0]), assembly_path, lefts,
           rights, cpu_cap, reference, transrate_dir)
    name = 'transrate_' + assembly_name
    out, err = GEN_LOGS(name)
    temp_task = Task(command=cmd, dependencies=[], targets=trgs, name=name, cpu=cpu_cap, stdout=out, stderr=err)
    deps = transrate_dep_generator(temp_task, orig_lefts, orig_rights, orig_singles, reference, assembly_path, cpu_cap, transrate_dir, tasks)
    temp_task.dependencies = [deps]
    return temp_task


def assembly_stats_task(out_dir,assembly,tasks):
    ''' Defines assembly_stats task. Uses PATH_DIR, PATH_SCRIPTS, NAME_ASSEMBLY.
        Params :
            tasks - a list of tasks that this task is dependant on (trinity_task)
    '''
    trgs = ['{0!s}/assembly_stats.json'.format(out_dir)]
    cmd = 'python {0!s}/assembly_stats.py {1!s} > {2!s}'.format(PATH_SCRIPTS,assembly,trgs[0])
    name = 'assembly_stats'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def build_bowtie_task(path_assembly, assembly_name, out_dir, tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.1.bt2'.format(out_dir,assembly_name)]  
    cmd = 'bowtie2-build --offrate 1 -f {1!s} {2!s}/{3!s}'.format(
            PATH_BOWTIE2, path_assembly, out_dir, assembly_name) 
    name = 'build_bowtie'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def bowtie2_unpaired_task(bowtie2_index,out_dir,fastq,out_name,opt,cpu_cap,tasks):
    '''
    '''
    opts = ['-a -t --end-to-end', '-t --local']
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir,out_name)]
    cmd = ('bowtie2 {1!s} -L {2!s} -N 1 --threads {3!s} -x {4!s} -U '
            '{5!s} | samtools view -Sb - > {6!s} ').format(PATH_BOWTIE2,
            opts[opt],22,cpu_cap,bowtie2_index,fastq,trgs[0])
    name = 'bowtie2_'+out_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def bowtie2_task(bowtie2_index,out_dir,fastq1,fastq2,out_name,opt,cpu_cap,tasks):
    '''    
    '''
    opts = ['-a -t --end-to-end', '-t --local']
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir,out_name)]
    cmd = ('bowtie2 {1!s} -L {2!s} -N 1 --maxins 800 --threads {3!s} -x {4!s} -1 '
            '{5!s} -2 {6!s} | samtools view -Sb - > {7!s} ').format(PATH_BOWTIE2,
            opts[opt],22,cpu_cap,bowtie2_index,fastq1,fastq2,trgs[0])
    name = 'bowtie2_'+out_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def express_task(assembly_path,out_dir,out_name,bam_input,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.xprs'.format(out_dir,out_name)]
    cmd = ('mkdir {1!s}/{2!s}; {0!s} --output-dir {1!s}/{2!s} {3!s} {4!s}; mv '
            '{1!s}/{2!s}/results.xprs {5!s}; rm -rf {1!s}/{2!s};').format(
            PATH_EXPRESS,out_dir,out_name,assembly_path,bam_input,trgs[0])
    name = 'express_'+out_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def counts_to_table_task(gene_trans_map,out_dir,express_files,out_name,flag,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.countsTable'.format(out_dir,out_name)]
    count_str = ' '.join(['--counts {0!s}'.format(f) for f in express_files])
    cmd = ('python {0!s}/counts_to_table2.py --out {1!s} --inDir {2!s} '
            '--outDir {2!s} {3!s} {4!s} --geneTransMap {5!s}').format( 
            PATH_SCRIPTS,out_name,out_dir,flag,count_str,gene_trans_map)
    name = 'counts_to_table_'+flag
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def sam_sort_task(out_dir,bam_file,out_name,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.bam'.format(out_dir,out_name)]
    cmd = 'samtools sort {0!s} {1!s}/{2!s}'.format(bam_file,out_dir,out_name)
    name = 'sam_sort_'+out_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def assembly_to_bed_task(path_assembly,assembly_name,out_dir, tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.bed'.format(out_dir,assembly_name)]
    cmd = 'python {0!s}/fasta_to_bed_count_length.py {1!s} {2!s}'.format(
            PATH_SCRIPTS,path_assembly,trgs[0])
    name = 'fasta_to_bed'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def intersect_bed_task(out_dir,bam_file,bed_reference,output_name,tasks):
    '''
    '''
    trgs = ['{0!s}/{1!s}.bed'.format(out_dir,output_name)]
    cmd = '{0!s} intersect -abam {1!s} -b {2!s} -wb -bed > {3!s}'.format(PATH_BEDTOOLS,bam_file,bed_reference,trgs[0])
    name = 'intersect_bed_'+output_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def deseq2_task(out_dir,counts_to_table_results,sample_info,basename,model,tasks):
    '''
    '''
    pseudo_model_temp = ''.join([c if c!=' ' else '_' for c in model])
    trgs = ['{0!s}/deseq2_{1!s}_{2!s}/'.format(out_dir,basename,pseudo_model_temp)]
    cmd = 'Rscript {5!s}/deseq2.r --args {0!s} {1!s} {2!s} {3!s} {4!s}'.format(
            counts_to_table_results,sample_info,out_dir,basename,model,PATH_SCRIPTS)
    name = 'de_'+basename
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def build_salmon_task(path_assembly,assembly_name,out_dir,cpu_cap,tasks):
    trgs = ['{0!s}/{1!s}_salmon'.format(out_dir, assembly_name)] 
    cmd = '{0!s} index --transcripts {1!s} --index {2!s}/{3!s}_salmon --threads {4!s} --type quasi'.format(tool_path_check(TOOLS_DICT['salmon'].full_exe[0]),path_assembly, out_dir, assembly_name, cpu_cap)
    name = 'build_salmon_' + assembly_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err, cpu=cpu_cap)


def salmon_gene_map_task(out_dir,assembly_name,gene_trans_map,tasks):
    ''' salmon requires gene_trans_map in reverse column order (transcript \t gene \n)'''
    trgs = ['{0!s}/{1!s}.trans_gene_map'.format(out_dir,assembly_name)] 
    cmd = 'join -t, -o 1.2,1.1 {0!s} {0!s} > {1!s}'.format(gene_trans_map, trgs[0]) 
    name = 'trans_gene_map'
    name = 'salmon_gene_map_task_' + assembly_name
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def salmon_task(index,left,right,out_name,gene_map,out_dir,cpu_cap,tasks):
    trgs = ['{0!s}/{1!s}/quant.sf'.format(out_dir,out_name)]
    cmd = '{0!s} quant -i {1!s} -l IU -1 {2!s} -2 {3!s} -o {4!s}/{5!s} --geneMap {6!s} -p {7!s} --extraSensitive; cp ' \
        '{4!s}/{5!s}/quant.sf {4!s}/{5!s}_quant.sf; cp {4!s}/{5!s}/quant.genes.sf {4!s}/{5!s}_quant.genes.sf'.format(
	tool_path_check(TOOLS_DICT['salmon'].full_exe[0]),index,left,right,out_dir,out_name,gene_map,cpu_cap)
    name = 'salmon_' + os.path.basename(index)
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def salmon_unpaired_task(index,unpaired,out_name,gene_map,out_dir,cpu_cap,tasks):
    trgs = []
    cmd = '{0!s} quant -i {1!s} -l U -r {2!s} -o {3!s}/{4!s} --geneMap {5!s} -p {6!s} --extraSensitive'.format(
            tool_path_check(TOOLS_DICT['salmon'].full_exe[0]),index,unpaired,out_dir,out_name,gene_map,cpu_cap)
    name = 'salmon_unpaired_' + os.path.basename(index)
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err,cpu=cpu_cap)


def build_kallisto_task(assembly_path, assembly_name,out_dir,tasks):
    trgs = []
    cmd = '{0!s} index -i {1!s}/{2!s}_kallisto {3!s}'.format(
            PATH_KALLISTO,out_dir,assembly_name,assembly_path)
    name = 'build_kallisto'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def kallisto_task(index,out_dir,out_name,left,right,tasks):
    trgs = []
    cmd = '{0!s} quant -i {1!s} -o {2!s}/{3!s} {4!s} {5!s}'.format(
            PATH_KALLISTO,index,out_dir,out_name,left,right)
    name = 'kallisto'
    out,err = GEN_LOGS(name)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def build_blast_task(fasta,out_path,dbtype,tasks,log_flag=True):
    trgs = []
    #title doesn't seem to change the out name .. it's still xx.gz.psq, etc? CHECK.
    title = os.path.basename(out_path).split('.')[0]
    cmd = 'gunzip -c {0!s} | makeblastdb -in - -dbtype {2!s} -title {3!s} -out {1!s}'.format(fasta,out_path,dbtype,title)
    name = 'build_blast_'+title
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def build_diamond_task(fasta,out_path,tasks,log_flag=True):
    title = os.path.basename(out_path)
    trgs = ['{0!s}'.format(out_path + '.dmnd')] 
    cmd = '{0!s} makedb --in {1!s} --db {2!s}'.format(TOOLS_DICT['diamond'].full_exe[0], fasta, out_path)
    name = 'build_diamond_'+ title
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)
 

def split_mito_task(blast_mt,tasks):
    trgs = ['{0!s}/mtDNA_contigs.fasta'.format(GEN_PATH_ASSEMBLY_FILES()),'{0!s}/no_mtDNA_contigs.fasta'.format(GEN_PATH_ASSEMBLY_FILES())]
    cmd = '{0!s}/split_fasta.py {1!s} {3!s} {2!s}/mtDNA_contigs.fasta {2!s}/no_mtDNA_contigs.fasta'.format(PATH_SCRIPTS,GEN_PATH_ASSEMBLY(),GEN_PATH_ASSEMBLY_FILES(),blast_mt)
    name = 'split_mito'
    out, err = GEN_LOGS(name) 
    return Task(command=cmd,dependencies=tasks,targets=trgs,name=name,stdout=out,stderr=err)


def pfam_build_task(source, tasks, log_flag=True):
    trgs = [PATH_PFAM_DATABASE+'.h3f']
    cmd = 'cd {0!s} ; {1!s} -f {2!s};'.format(PATH_DATABASES, tool_path_check(TOOLS_DICT['hmmer'].full_exe[1]), source)
    name = 'hmmpress'
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)


def db2stitle_task(db, tasks, log_flag=True):
    base_db = os.path.basename(db)
    trgs = ['{0!s}/{1!s}.stitle'.format(PATH_DATABASES, base_db)]
    cmd = 'python {0!s}/fastaID2names.py --fasta {1!s} > {2!s}'.format(
           PATH_SCRIPTS, db, trgs[0])
    name = 'db2stitle_'+base_db
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)

def manage_db_task(fresh, nr_flag, uniref90_flag, busco_flags, blastplus_flag, cpu_cap, tasks, log_flag=True):
    trgs = [PATH_DATABASES]
    cmd = 'python {0!s}/manage_database.py'.format(PATH_SCRIPTS)
    if(fresh):
        cmd +=' --hard'
    if(nr_flag):
        cmd+=' --nr'
    if(uniref90_flag):
        cmd+=' --uniref90'
    if(len(busco_flags) != 0):
        cmd+=' --buscos '
        cmd+=','.join(busco_flags)
    if blastplus_flag:
        cmd+= ' --buildBlastPlus'
    name = 'db_manage'
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap)

def manage_tools_task(install, fresh, cpu_cap, tool_list, tasks, log_flag=True):
    trgs = [PATH_TOOLS]
    cmd = 'python {0!s}/manage_tools.py'.format(PATH_SCRIPTS)
    if(install):
        cmd+=' --install'
    if(fresh):
        cmd +=' --hard'
    cmd+= ' --tool ' + ' --tool '.join(tool_list)
    name = 'tools_manage'
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap)

def filter_task(assembly_path, assembly_name, out_dir, quant_file_list, tpm_threshold, tpm_column_index, tasks, log_flag=True):
    # TPM column index: transrate uses older salmon; use index =2. Newer salmon: index=3
    trgs = ['{0!s}/{1!s}_{2!s}tpm.fasta'.format(out_dir,assembly_name,tpm_threshold)]
    quants = ''.join(' --quant_files '+ x for x in quant_file_list) 
    cmd = 'python {0!s}/filter_contigs_by_tpm.py --assembly {1!s} --tpm {2!s} {3!s} --out {4!s} --tpm_column_index {5!s}'.format(PATH_SCRIPTS, assembly_path,tpm_threshold, quants, trgs[0], tpm_column_index)
    name = 'filt_{0!s}_{1!s}tpm'.format(assembly_name, tpm_threshold)
    out, err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)

