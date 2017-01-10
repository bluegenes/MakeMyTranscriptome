'''
'''
from tasks_v2 import Task, Supervisor
import os
from external_tools import TOOLS_DICT
from functions_general import gen_logs, tool_path_check, make_dir_task, cp_task
import mmt_defaults as statics

# opc is the output_path_class object


def fastqc_task(opc, out_dir, fq_files, output_name, cpu_cap, tasks):
    '''    Defines task for running fastqc. Uses GEN_PATH_DIR(), PATH_FASTQC,
        Params :
            fq_files - list of fastq files to run fastqc on
            tasks - a list of tasks that this task is dependent on.
    '''
    cpu_param = min(len(fq_files), cpu_cap)
    outDir = '{0!s}/fastqc_{1!s}'.format(out_dir, output_name)
    trgs = [os.path.join(outDir, os.path.basename(x).rsplit('.f')[0] + '_fastqc.zip') for x in fq_files]
    cmd = 'mkdir {2!s}; {0!s} --extract --outdir {2!s} --threads {3!s} {1!s}'.format(
           TOOLS_DICT['fastqc'].full_exe[0],' '.join(fq_files), outDir, cpu_param)
    name = 'fastqc_'+output_name
    out, err = gen_logs(opc.path_logs, name)
    fast_task = Task(command=cmd, dependencies=[mkdir_task], targets=trgs, name=name, stdout=out, stderr=err)
    super_name = "super_" + name
    return Supervisor([mkdir_task, fast_task], name=super_name, dependencies=tasks)


def prinseq_unpaired_task(opc, out_dir, input1, basename, opts, tasks):
    trgs = ['{0!s}/{1!s}_{2!s}'.format(out_dir, basename, os.path.basename(input1))]
    pseudo_trgs =['{0!s}/{1!s}.fastq'.format(out_dir, basename)]
    cmd = ('perl {0!s} -fastq {1!s} --out_format 3 --out_good {2!s}/{3!s} --out_bad null '
           '--trim_qual_left 20 --trim_qual_right 20 --trim_qual_type min --min_len 35 '
           '--trim_tail_left 8 --trim_tail_right 8 {4!s} -log'
           ).format(tool_path_check(TOOLS_DICT['prinseq'].full_exe[0]), input1, out_dir,
                    basename, opts, trgs[0])
    name = basename
    out, err = gen_logs(opc.path_logs, name)
    prinseq = Task(command=cmd, dependencies=[], name=name, stdout=out, stderr=err, targets=pseudo_trgs)
    cp = cp_task(pseudo_trgs[0], trgs[0], [prinseq])
    super_name = "super_" + name
    return Supervisor([prinseq, cp], dependencies=tasks, name=super_name)


def prinseq_task(opc, out_dir, input_1, input_2, basename, opts, tasks):
    '''    Defines prinseq task. Uses GEN_PATH_DIR(), PATH_PRINSEQ
        Params :
            input_1 - a list of 1/left fastq files
            input_2 - a list of 2/right fastq files
            basename - the basename for all output files
            opts - optional params for trinity task.
            tasks = the tasks that this task is dependent on
    '''
    trgs = ['{0!s}/{1!s}_1_{2!s}'.format(out_dir, basename, os.path.basename(input_1)),
            '{0!s}/{1!s}_2_{2!s}'.format(out_dir, basename, os.path.basename(input_2))]
    pseudo_trgs = ['{0!s}/{1!s}_{2!s}.fastq'.format(out_dir, basename, x) for x in range(1, 3)]
    cmd = ('perl {0!s} -fastq {1!s} -fastq2 {2!s} --out_format 3 --out_good {3!s}/{4!s} '
           '--out_bad null --trim_qual_left 20 --trim_qual_right 20 --trim_qual_type min '
           '--min_len 55 --trim_tail_left 8 --trim_tail_right 8 {5!s} -log'
           ).format(tool_path_check(TOOLS_DICT['prinseq'].full_exe[0]), input_1,
           input_2, out_dir, basename, opts, pseudo_trgs[0], trgs[0], pseudo_trgs[1], trgs[1])
    name = basename
    out, err = gen_logs(opc.path_logs, name)
    prinseq = Task(command=cmd, dependencies=[], name=name, stdout=out, stderr=err, targets=pseudo_trgs)
    cp1 = cp_task(pseudo_trgs[0], trgs[0], [prinseq])
    cp2 = cp_task(pseudo_trgs[1], trgs[1], [prinseq])
    super_name = "super_" + name
    return Supervisor(tasks=[prinseq, cp1, cp2], dependencies=[tasks], name=super_name)


def trimmomatic_unpaired_task(opc, out_dir,input1, cpu_cap, basename, tasks):
    trgs = ['{0!s}/{1!s}_{2!s}'.format(out_dir, basename, os.path.basename(input1))]
    cmd = ('java -jar {0!s} SE -threads {3!s} {1!s} {2!s} ILLUMINACLIP:'
           '{4!s}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35'
           ).format(tool_path_check(TOOLS_DICT['trimmomatic'].full_exe[0]),
           input1, trgs[0],cpu_cap, TOOLS_DICT['trimmomatic'].full_exe[2])  # PATH_TRIMMOMATIC_ADAPTERS_SINGLE)
    name = basename
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs, cpu=cpu_cap) 


def trimmomatic_task(opc, out_dir, left, right, cpu_cap, basename, tasks):
    base_str = '{0!s}/{1!s}'.format(out_dir, basename)
    trgs = [base_str+'_1_' + os.path.basename(left),
            base_str+'_2_' + os.path.basename(right)]
    orphans = [base_str+'_1s_' + os.path.basename(left),
               base_str+'_2s_' + os.path.basename(right)]
    cmd = ('java -jar {0!s} PE -threads {3!s} {1!s} {2!s} {5!s} {4!s} {7!s} '
           '{6!s} ILLUMINACLIP:{8!s}:2:30:10 LEADING:3 TRAILING:3 '
           'SLIDINGWINDOW:4:15 MINLEN:35').format(
           tool_path_check(TOOLS_DICT['trimmomatic'].full_exe[0]), left, right,
           cpu_cap, orphans[0], trgs[0], orphans[1], trgs[1],
           TOOLS_DICT['trimmomatic'].full_exe[1])
    name = basename
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs, cpu=cpu_cap) 


def rcorrector_task(opc, out_dir, left, right, cpu_cap, basename, tasks):
    trgs = ['{0!s}/{1!s}.cor.fq'.format(out_dir, os.path.basename(left)),
            '{0!s}/{1!s}.corr.fq'.format(out_dir, os.path.basename(right))]
    cmd = 'perl {0!s} -1 {1!s} -2 {2!s} -t {3!s}'.format(
          tool_path_check(TOOLS_DICT['rcorrector'].full_exe[0]), left, right, cpu_cap)
    name = 'Rcorrector_' + basename
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs, cpu=cpu_cap)


def remove_dups_task(opc, out_dir, left, right, out_base, tasks):
    '''    Defines rmdup_fastq_paired tasks. Uses PATH_SCRIPTS, GEN_PATH_DIR()
        Params : 
            left : a set of left/1 files to have duplicates removed from
            right : a set of right/2 files to have duplicates removed from
            out_base : Basename for output files
            tasks : A set of tasks that this task is dependent on.
    '''
    trgs = ['{0!s}/{1!s}_1.fastq'.format(out_dir, out_base),
            '{0!s}/{1!s}_2.fastq'.format(out_dir, out_base)]
    cmd = ('python {0!s}/rmdup_fastq_paired.py --left {1!s} --right {2!s} '
           '--left_target {3!s} --right_target {4!s}').format(
           statics.PATH_UTIL, ','.join(left), ','.join(right), trgs[0], trgs[1])
    name = 'remove_dups'
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs)


def basic_cat_task(sources, target, tasks):
    trgs = [target]
    cmd = 'cat {0!s}'.format(' '.join(sources))
    out = target
    err = os.devnull
    name = "cat_" + os.path.basename(target)
    return Task(command=cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs)


def cat_task(opc, out_dir, left, right, basename, tasks):
    trgs = ['{0!s}/{1!s}_1.fastq'.format(out_dir, basename),
            '{0!s}/{1!s}_2.fastq'.format(out_dir, basename)]
    left_cat = basic_cat_task(left, trgs[0], [])
    right_cat = basic_cat_task(right, trgs[1], [])
    name = 'cat_' + basename
    return Supervisor(tasks=[left_cat, right_cat], name=name, dependencies=[tasks])


def subset_task(opc, out_dir, fastq1, fastq2, out_base, num, seed, tasks):
    trgs = ['{0!s}/{1!s}_1.fastq'.format(out_dir, out_base),
            '{0!s}/{1!s}_2.fastq'.format(out_dir, out_base)]
    cmd = 'python {0!s}/random_subset.py -1 {1!s} -2 {2!s} -n 100 -s {3!s} -t1 {4!s} -t2 {5!s}'.format(
            statics.PATH_UTIL, ','.join(fastq1), ','.join(fastq2), num, trgs[0], trgs[1])
    if(seed is not None):
        cmd += ' --seed '+seed
    name = 'subset_reads'
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs)

#def seqtk_subset_task(out_dir,left, right, num_seqs, seed,tasks):#out_dir, fastq1, fastq2, out_base, num, seed, tasks):
#     form = lambda s, i : s.format(out_dir, os.path.basename(i), num_seqs)
#     trgs = [form('{0!s}/{1!s}_sub{2!s}.fq', left),
#             form('{0!s}/{1!s}_sub{2!s}.fq', right)]
#    cmd = '{0!s} sample -s {1!s} {2!s} > {3!s}; {0!s} sample -s {1!s} {4!s} > {5!s};'.format(
#           fg.tool_path_check(TOOLS_DICT['seqtk'].full_exe[0],seed,left,num_seqs,trgs[0],right,trgs[1])
#    name = 'seqtk_' + basename
#    out, err = fg.GEN_LOGS(name)
#    return Task(command=cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs, cpu=cpu_cap)

#def seqtk_truncate_task(out_dir,left,right,cpu_cap,basename,tasks):
#    out, err = fg.GEN_LOGS(name)
#    return Task(command=cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=trgs, cpu=cpu_cap)


def basic_truncate_task(opc, source, target, length, tasks):
    cmd = 'python {0!s}/truncate_fastq.py {1!s} --length {2!s} --target {3!s}'.format(
          statics.PATH_UTIL, source, length, target)
    name = "truncate_" + os.path.basename(source)
    out, err = gen_logs(opc.path_logs, name)
    return Task(command=cmd, dependencies=tasks, name=name, stdout=out, stderr=err, targets=[target])


def truncate_task(opc, out_dir, left, right, length, tasks):
    trgs = ['{0!s}/truncated_1.fastq'.format(out_dir, left),
            '{0!s}/truncated_2.fastq'.format(out_dir, right)]
    trunc1 = basic_truncate_task(opc, left, trgs[0], length, [])
    trunc2 = basic_truncate_task(opc, right, trgs[1], length, [])
    name = 'truncate_reads'
    return Supervisor(tasks=[trunc1, trunc2], name=name, dependencies=tasks)


def trinity_task(opc, path_assembly, out_dir, fastq, fastq2, unpaired, cpu_cap_trin, cpu_cap_bfly, mem_trin, mem_bfly, normalize_flag, tasks):
    '''    Defines the trinity task. Uses GEN_PATH_DIR(), PATH_TRINITY, NAME_ASSEMBLY
        Params :
            left - a 1/left fastq files
            right - a 2/right fastq files
            cpu_cap - number of threads used by trinity
            tasks - a list of tasks that this task is dependent on
    '''
    normalize_flag = '--normalize_reads' if(normalize_flag) else ''
    input_str = ''
    if(unpaired != [] and fastq == []):
        input_str += '--single ' + ','.join(unpaired)
    if(fastq != []):
        input_str += '--left ' + ','.join(fastq + unpaired)
        input_str += ' --right ' + ','.join(fastq2)
    trgs = [path_assembly]
    pseudo_trgs = ['{0!s}/trinity/Trinity.fasta'.format(out_dir)]
    cmd = ('{0!s} --seqType fq {1!s} --CPU {2!s} --max_memory {3!s}G --bflyCalculateCPU {4!s} '
           '--output {6!s}/trinity'
           ).format(tool_path_check(TOOLS_DICT['trinity'].full_exe[0]), input_str, cpu_cap_trin,
                    mem_trin, normalize_flag, mem_bfly, out_dir, trgs[0])
    name = 'trinity_assembly'
    out, err = gen_logs(opc.path_logs, name)
    cpu_cap = max(cpu_cap_trin, cpu_cap_bfly)
    trinity = Task(command=cmd, dependencies=[], targets=pseudo_trgs, name=name, cpu=cpu_cap, stdout=out, stderr=err)
    cp = cp_task(pseudo_trgs[0], trgs[0], [trinity])
    super_name = 'super_' + name
    return Supervisor(tasks=[trinity, cp], dependencies=tasks, name=super_name)


def rnaspades_task(opc, path_assembly, out_dir, left, right, unpaired, cpu_cap, tasks):
    virtual_target = '{0!s}/rna_spades_out_dir'.format(out_dir)
    trgs = [path_assembly]
    pseudo_trgs = ['{0!s}/contigs.fasta'.format(virtual_target)]
    input_strings = []
    if(left != []):
        input_strings.append('-1 '+left[0])
        input_strings.append('-2 '+right[0])
    if(unpaired != []):
        input_strings.append('-s '+unpaired[0])
    cmd = '{0!s} {1!s} --threads {2!s} -o {3!s}'.format(
            tool_path_check(TOOLS_DICT['rnaspades'].full_exe[0]), ' '.join(input_strings),
            cpu_cap, virtual_target, trgs[0])
    name = 'rnaSPAdes_assembly'
    out, err = gen_logs(opc.path_logs, name)
    rnaspades = Task(command=cmd, dependencies=[], targets=trgs, name=name, stdout=out, stderr=err, cpu=cpu_cap)
    cp = cp_task(pseudo_trgs[0], trgs[0], [rnaspades])
    super_name = 'super_' + name
    return Supervisor(tasks=[rnaspades, cp], dependencies=task, name=super_name)
