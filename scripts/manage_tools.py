from time import strftime
import tarfile
from task_functions_v2 import PATH_ROOT
import os
import json
import argparse
import sys
import functools
import subprocess
import platform
from tasks_v2 import Supervisor
if(sys.version[0] == '3'):
    from urllib.request import urlretrieve, ContentTooShortError
    import shutil
    which = shutil.which
else:
    from urllib import urlretrieve, ContentTooShortError
    which = which_python2


trinity_linux_url = 'https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz'
trinity_linux_target = os.path.join(PATH_TOOLS, 'trinityrnaseq-2.1.1')
trinity_exe = 'Trinity'

trimmomatic_url = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip'
trimmomatic_trinity_target = os.path.join(which('Trinity'), 'trinity-plugins/Trimmomatic')
trimmomatic_url_target = os.path.join(PATH_TOOLS, 'Trimmomatic-0.35')
trimmomatic_exe = 'Trimmomatic-0.35.jar'

prinseq_url = 'http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz'
prinseq_target = os.path.join(PATH_TOOLS, 'prinseq-lite-0.20.4')
prinseq_exe = 'prinseq-lite.pl'

transdecoder_url = 'https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz'
transdecoder_target = os.path.join(PATH_TOOLS, 'TransDecoder-2.0.1')
transdecoder_exe1 = 'TransDecoder.Predict'
transdecoder_exe2 = 'TransDecoder.LongOrfs'

transrate_linux_url = 'https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz'
transrate_linux_target = os.path.join(PATH_TOOLS, 'transrate-1.0.1-linux-x86_64')
transrate_exe = 'transrate'

busco_url = 'http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz'
busco_target = os.path.join(PATH_TOOLS, 'BUSCO_v1.1b1')
busco_exe = 'BUSCO_v1.1b1.py'

hmmer_linux_url = 'http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz'
hmmer_linux_target = os.path.join(PATH_TOOLS, 'hmmer-3.1b2-linux-intel-x86_64')
hmmer_exe1 = 'hmmscan' 
hmmer_exe2 = 'hmmpress'

diamond_linux_url = 'http://github.com/bbuchfink/diamond/releases/download/v0.7.10/diamond-linux64.tar.gz'
diamond_linux_target = os.path.join(PATH_TOOLS, 'diamond')
diamond_exe = 'diamond'

salmon_linux_url = 'https://github.com/COMBINE-lab/salmon/releases/download/v0.5.1/SalmonBeta-0.5.1_DebianSqueeze.tar.gz'
salmon_linux_target = os.path,join(PATH_TOOLS, 'SalmonBeta-0.5.1_DebianSqueeze')
salmon_exe = 'salmon'

#	trinity_version_check = 'Trinity --version'
#	prinseq_version_check = 'perl prinseq-lite.pl -version'
#	transrateversion_check = 'transrate --version'

def which_python2(program):
    #from: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python/377028#377028
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

tool_supervisor_log = '{0!s}/.tool_supervisor_log'.format(PATH_TOOLS)

def check_tools():
    if(not os.path.isdir(PATH_TOOLS)):
        os.mkdir(PATH_TOOLS)
    sys.path.append(PATH_TOOLS) 
    if(which(trinity_exe):
	toolsD['trinity'] = True
    if(which(prinseq_exe)):
	toolsD['prinseq'] = True
    if(which(transrate_exe):
        toolsD['transrate'] = True
    if(which(busco_exe)):
        toolsD['busco'] = True
    if(which(transdecoder_exe1) and which(transdecoder_exe2)):
	toolsD['transdecoder'] = True
    if(which(diamond_exe)):
        toolsD['diamond'] = True
    if(which(hmmer_exe1) and which(hmmer_exe2):
        toolsD['hmmer'] = True
    if(which(salmon_exe)):
        toolsD['salmon'] = True
    if(not os.path.isfile(os.path.join(PATH_TOOLS, '.tools_log'))):
        write_log({})


def print_install_instructions(toolsDt, syst):
    instructions = "some instructions will go here"
    # global variable = links, installation info for each program? or pass in a dictionary
    print(instructions)

def download_external_tools(log_table, file_check=True):
    partial_get = lambda a, b, c : get(log_table, a, b ,c, file_check)
#    partial_get('', url_go_pathway, go_pathway_target)
    #if(uniref90_flag):
    #    partial_get('', url_uniref90, uniref90_target)
    #partial_get('gz', url_id_mapping, id_mapping_target)
    #if(busco_flags['metazoa']):
    #    partial_get('tar', url_busco_metazoa, busco_metazoa_target)
    return log_table


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


#def version_check(tool_name, version_command):
#    print('checking for ' + tool_name)
#    try:
#        subprocess.call([tool_name, version_command])
#    except OSError as e:
#        if e.errno == os.errno.ENOENT:
        # handle file not found error.
#        else:
        # Something else went wrong while trying to run `wget`
#            raise


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
    log = open(os.path.join(PATH_TOOLS, 'tool_log'))
    log_table = json.load(log)
    log.close()
    return log_table


def write_log(log_table):
    log = open(os.path.join(PATH_DATABASES, 'database_log'), 'w')
    json.dump(log_table, log, sort_keys=True, indent=4)
    log.close()


def main(install=False, reinstall=False, trimmomatic = False, cpu=4):
    toolsD = {'trinity' = False, 'trimmomatic' = False, 'prinseq' = False, 'transdecoder' = False, 'transrate' = False, 'busco' = False, 'diamond' = False, 'hmmer' = False, 'salmon' = False}
    toolsD = check_tools(toolsD) #return altered tools dictionary
    #optionalToolsDt = {} #add non-required tool checks: bowtie2, express, bedtools, etc
    log_table = read_log()
#    log_table = download_tools(log_table, toolsD)
    tasks = []
    if(install and platform.system().lower() == 'linux'):
        if(not toolsD['trinity']):
	    download_external_tool(trinity_linux_url, trinity_linux_target)
            install_trinity = install_trinity_task(trinity_linux_target, trinity_exe)
	    tasks.append(install_trinity)
	if(not toolsD['prinseq']):
	    download_external_tool(prinseq_url, prinseq_target)
	    install_prinseq = install_prinseq_task(prinseq_target, prinseq_exe)
            tasks.append(install_prinseq)
	if(not toolsD['transdecoder']):
	    download_external_tool(transdecoder_url, transdecoder_target)
	    install_transdecoder = install_transdecoder_task(transdecoder_target, transdecoder_exe)
	    tasks.append(install_transdecoder)
	if(not toolsD['transrate']):
	    download_external_tool(transrate_linux_url, transrate_linux_target)
	    install_transrate = install_transrate_task(transrate_linux_target, transrate_exe)
	    tasks.append(install_transrate)
	if(not toolsD['busco']):
	    download_external_tool(busco_url, busco_target)
	    install_busco = install_busco_task(busco_target, busco_exe)
	    tasks.append(install_busco)
	if(not toolsD['diamond']):
	    download_external_tool(diamond_linux_url, diamond_linux_target)
	    install_diamond = install_diamond_task(diamond_linux_target, diamond_exe)
	    tasks.append(install_diamond)
	if(not toolsD['hmmer']):
	    download_external_tool(hmmer_linux_url, hmmer_linux_target)
	    install_hmmer = install_hmmer_task(hmmer_linux_target, hmmer_exe1, hmmer_exe2)
	    tasks.append(install_hmmer)
	if(not toolsD['salmon']):
	    download_external_tool(salmon_linux_url, salmon_linux_target)
	    install_salmon = install_salmon_task(salmon_linux_target, salmon_exe)
	    tasks.append(install_salmon)
#        if(not toolsD['trimmomatic']):
#	    download_external_tool(prinseq_url, trimmomatic_url_target)
#            install_trimmomatic = install_trimmomatic_task(trimmomatic_url_target, trimmomatic_exe)
#	    tasks.append(install_trimmomatic)
    else:
        print_install_instructions(toolsD)
    run_tasks(tasks, cpu)
    write_log(log_table)


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument('--install', action='store_true', default=False)
    parser.add_argument('--hard', action='store_true', default=False)
    parser.add_argument('--cpu', type=int)
    args = parser.parse_args()
    main(args.install, not args.hard, args.cpu)

