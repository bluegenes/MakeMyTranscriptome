from time import strftime
import tarfile
from task_functions_v2 import PATH_ROOT
import os
import json
import argparse
import sys
import functools
import subprocess
from tasks_v2 import Supervisor
if(sys.version[0] == '3'):
    from urllib.request import urlretrieve, ContentTooShortError
else:
    from urllib import urlretrieve, ContentTooShortError

if not os.path.exists(PATH_TOOLS):
    os.makedirs(PATH_TOOLS)

tool_supervisor_log = '{0!s}/.tool_supervisor_log'.format(PATH_TOOLS)

# need to change all targets --> what they need to be
#ASSEMBLY TOOLS#
url_trinity = 'https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz'
trinity_target = os.path.join(PATH_TOOLS, 'trinity')
# need to cd into Trinity dir; 'make' 
#trinity_version_check = 'Trinity --version'
# both trimmomatic and transdecoder are included as part of trinity --can just link to those versions! 

#url_trimmomatic = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip'
#trimmomatic_target = os.path.join(PATH_TOOLS, 'trimmommatic')

trimmomatic_version_check = 
#url_prinseq = 'http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz'
#prinseq_target = os.path.join(PATH_TOOLS, 'prinseq')
prinseq_version_check = 'prinseq -version'
#QUALITY TOOLS#
#url_transrate_linux = 'https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz'
#url_transrate_mac = 'https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-osx.tar.gz'
#transrate_target = os.path.join(tools_folder,'transrate')
transrate_version_check = 'transrate --version'
#url_busco = 'http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz'
#busco_target = os.path.join(tools_folder,'busco')

#ANNOTATION TOOLS#
#url_transdecoder = 'https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz'
#transdecoder_target = os.path.join(tools_folder, 'transdecoder')
#url_diamond_linux = 'http://github.com/bbuchfink/diamond/releases/download/v0.7.10/diamond-linux64.tar.gz'
#url_diamond_source = 'http://github.com/bbuchfink/diamond/archive/v0.7.10.tar.gz' 
#diamond_target = os.path.join(tools_folder, 'diamond')
#url_hmmer = 'http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz'
#hmmer_target = os.path.join(tools_folder, 'hmmer')
#EXPRESSION TOOLS#
#url_salmon_linux = 'https://github.com/COMBINE-lab/salmon/releases/download/v0.5.1/SalmonBeta-0.5.1_DebianSqueeze.tar.gz'
#url_salmon_source = 'https://github.com/COMBINE-lab/salmon/archive/v0.5.1.tar.gz'
#salmon_target = os.path,join(tools_folder, 'salmon')




def run_tasks(tasks, cpu=4):
    for t in tasks:
        print(t.name)
	t.stdout = os.path.join(PATH_TOOLS, t.name, '.stdout')
        t.stderr = os.path.join(PATH_TOOLS, t.name,'.stderr')

    s = Supervisor(tasks=tasks, force_run=False, log=database_supervisor_log, cpu=cpu)
    s.run()
    for t in tasks:#if everything executes properly, rm the task logs
        if os.path.exists(t.stdout):
	    os.remove(t.stdout)
        if os.path.exists(t.stderr):
            os.remove(t.stderr)


def safe_retrieve(source, target):
    print('getting '+source)
    urlretrieve(source, target+'.temp')
    os.rename(target+'.temp', target)


def url_unzip(source, target):
    print('getting '+source)
    urlretrieve(source, target+'.gz')
    f = gzip.open(target+'.gz', 'rb')
    g = open(target, 'wb')
    for line in f:
        g.write(line)
    f.close()
    g.close()
    os.remove(target+'.gz')


def tar_retrieve(source, target):
    print('getting '+source)
    urlretrieve(source, target+'.tar.gz')
    tfile = tarfile.open(target+'.tar.gz', 'r:gz')
    tfile.extractall(target)
    os.remove(target+'.tar.gz')

def check_path_for_tool(tool_name, command):
    print('checking for ' + tool_name)
    try:
        subprocess.call([tool_name, command])
    except OSError as e:
        if e.errno == os.errno.ENOENT:
        # handle file not found error.
        else:
        # Something else went wrong while trying to run `wget`
            raise


    command -v foo >/dev/null 2>&1 || { echo >&2 "I require foo but it's not installed.  Aborting."; exit 1; }
    name = 'check_path_' + tool_name
    out_err = GEN_LOGS(name) if(log_flag) else (None, None)
    return Task(command=cmd, dependencies=tasks, targets=trgs, name=name, stdout=out, stderr=err)

def get(log_table, flag, source, target, file_check=True):
    if(file_check and os.path.exists(target)):
        return
    try:
        if(flag == 'gz'):
            url_unzip(source, target)
        elif(flag == ''):
            safe_retrieve(source, target)
        elif(flag == 'tar'):
            tar_retrieve(source, target)
        else:
            print('Can\'t retrieve database.')
    except ContentTooShortError:
        print('failed to install {0!s}'.format(source))
    basename = os.path.basename(target)
    log_table[basename] = strftime('%b-%d-%Y')


def read_log():
    log = open(os.path.join(PATH_DATABASES, 'database_log'))
    log_table = json.load(log)
    log.close()
    return log_table


def write_log(log_table):
    log = open(os.path.join(PATH_DATABASES, 'database_log'), 'w')
    json.dump(log_table, log, sort_keys=True, indent=4)
    log.close()


def download_external_tools(log_table, file_check=True):
    partial_get = lambda a, b, c : get(log_table, a, b ,c, file_check)
    partial_get('', url_go_pathway, go_pathway_target)
    if(uniref90_flag):
        partial_get('', url_uniref90, uniref90_target)
    partial_get('gz', url_id_mapping, id_mapping_target)
    if(busco_flags['metazoa']):
        partial_get('tar', url_busco_metazoa, busco_metazoa_target)
    return log_table


def check_tools():
    if(not os.path.isdir(PATH_DATABASES)):
        os.mkdir(PATH_DATABASES)
    if(not os.path.isdir(swissprot_folder)):
        os.mkdir(swissprot_folder)
    if(not os.path.isdir(uniref90_folder)):
        os.mkdir(uniref90_folder)
    if(not os.path.isdir(nr_folder)):
        os.mkdir(nr_folder)
    if(not os.path.isdir(pfam_folder)):
        os.mkdir(pfam_folder)
    if(not os.path.isdir(busco_folder)):
        os.mkdir(busco_folder)
    if(not os.path.isfile(os.path.join(PATH_DATABASES, '.database_log'))):
        write_log({})


def main(nr_flag=False, uniref90_flag=False, file_check=True, busco_flags=busco_flags, blastplus=False, cpu=4):
    check_database_dir()
    log_table = read_log()
    log_table = download_databases(log_table, nr_flag, uniref90_flag, file_check, busco_flags)
    log_table = subset_dat(id_mapping_target, idmapping_keys, log_table)
    tasks = []
    if blastplus:
        swissprot_task = build_blast_task(sprot_target, sprot_target, 'prot', [], False)
        tasks.append(swissprot_task)
    swissprot_diamond = build_diamond_task(sprot_target, PATH_SWISS_PROT, [], False)
    tasks.append(swissprot_diamond)
    swissprot_table_task = db2stitle_task(sprot_target, [], False)
    tasks.append(swissprot_table_task)
    if(uniref90_flag and os.path.exists(uniref90_target)):
        uniref90_diamond = build_diamond_task(uniref90_target, PATH_UNIREF90, [], False)
        tasks.append(uniref90_diamond)
        uniref90_table_task = db2stitle_task(uniref90_target, [], False)
        tasks.append(uniref90_table_task)
        if blastplus:
  	    uniref90_task = build_blast_task(uniref90_target,  PATH_UNIREF90, 'prot', [], False)
            tasks.append(uniref90_task)
    if(nr_flag and os.path.exists(nr_target)):
        nr_diamond = build_diamond_task(nr_target, PATH_NR, [], False)
        tasks.append(nr_diamond)
	nr_table_task = db2stitle_task(nr_target, [], False)
        tasks.append(nr_table_task)
        if blastplus:
	    nr_task = build_blast_task(nr_target, PATH_NR, 'prot', [], False)
	    tasks.append(nr_task)
    pfam_task = pfam_build_task(pfam_db_target, [], False)
    tasks.append(pfam_task)
    run_tasks(tasks, 4) # NEED TO FIX CPU HERE?
    write_log(log_table)


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument('--hard', action='store_true', default=False)
    parser.add_argument('--uniref90', action='store_true', default=False)
    parser.add_argument('--nr', action='store_true', default=False)
    parser.add_argument('--buscos', help='a comma seperated list of busco files that need to be downloaded')
    parser.add_argument('--cpu', type=int)
    parser.add_argument('--buildBlastPlus', action='store_true', default=False)
    args = parser.parse_args()
    if(args.buscos != None):
    	args.buscos = args.buscos.split(',')
    	for b in args.buscos:
            busco_flags[b] = True
    main(args.nr, args.uniref90, not args.hard, busco_flags, args.buildBlastPlus, args.cpu)
