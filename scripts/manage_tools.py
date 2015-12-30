from time import strftime
import tarfile
from task_functions_v2 import (PATH_ROOT, PATH_TOOLS,
install_trinity_task, install_trimmomatic_task, install_prinseq_task,
install_transdecoder_task, install_hmmer_task, install_salmon_task )
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
    from py2_which import which_python2
    which = which_python2


trinity_linux_url = 'https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz'
trinity_linux_target = os.path.join(PATH_TOOLS, 'trinityrnaseq-2.1.1')
#trinity_linux_target = os.path.join(PATH_TOOLS, '_trinity')
trinity_exe = 'Trinity'

trimmomatic_url = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip'
trimmomatic_trinity_target = os.path.join(which('Trinity'), 'trinity-plugins/Trimmomatic')
trimmomatic_url_target = os.path.join(PATH_TOOLS, 'Trimmomatic-0.35')
#trimmomatic_url_target = os.path.join(PATH_TOOLS, '_trimmomatic')
trimmomatic_exe = 'Trimmomatic-0.35.jar'

prinseq_url = 'http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz'
prinseq_target = os.path.join(PATH_TOOLS, 'prinseq-lite-0.20.4')
#prinseq_target = os.path.join(PATH_TOOLS, '_prinseq')
prinseq_exe = 'prinseq-lite.pl'

transdecoder_url = 'https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz'
transdecoder_target = os.path.join(PATH_TOOLS, 'TransDecoder-2.0.1')
#transdecoder_target = os.path.join(PATH_TOOLS, '_transdecoder')
transdecoder_exe1 = 'TransDecoder.Predict'
transdecoder_exe2 = 'TransDecoder.LongOrfs'

transrate_linux_url = 'https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz'
transrate_linux_target = os.path.join(PATH_TOOLS, 'transrate-1.0.1-linux-x86_64')
#transrate_linux_target = os.path.join(PATH_TOOLS, '_transrate')
transrate_exe = 'transrate'

busco_url = 'http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz'
busco_target = os.path.join(PATH_TOOLS, 'BUSCO_v1.1b1')
#busco_target = os.path.join(PATH_TOOLS, '_busco')
busco_exe = 'BUSCO_v1.1b1.py'

hmmer_linux_url = 'http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz'
hmmer_linux_target = os.path.join(PATH_TOOLS, 'hmmer-3.1b2-linux-intel-x86_64')
#hmmer_linux_target = os.path.join(PATH_TOOLS, '_hmmer3')
hmmer_exe1 = 'hmmscan' 
hmmer_exe2 = 'hmmpress'

diamond_linux_url = 'http://github.com/bbuchfink/diamond/releases/download/v0.7.10/diamond-linux64.tar.gz'
diamond_linux_target = os.path.join(PATH_TOOLS, '_diamond')
diamond_exe = 'diamond'

salmon_linux_url = 'https://github.com/COMBINE-lab/salmon/releases/download/v0.5.1/SalmonBeta-0.5.1_DebianSqueeze.tar.gz'
salmon_linux_target = os.path.join(PATH_TOOLS, 'SalmonBeta-0.5.1_DebianSqueeze')
#salmon_linux_target = os.path.join(PATH_TOOLS, '_salmon')
salmon_exe = 'salmon'

#	trinity_version_check = 'Trinity --version'
#	prinseq_version_check = 'perl prinseq-lite.pl -version'
#	transrateversion_check = 'transrate --version'

tool_supervisor_log = '{0!s}/.tool_supervisor_log'.format(PATH_TOOLS)

def check_tools(toolsD):
    if(which(trinity_exe)):
	toolsD['trinity'] = True
    if(which(prinseq_exe)):
	toolsD['prinseq'] = True
    if(which(transrate_exe)):
        toolsD['transrate'] = True
    if(which(busco_exe)):
        toolsD['busco'] = True
    if(which(transdecoder_exe1) and which(transdecoder_exe2)):
	toolsD['transdecoder'] = True
    if(which(diamond_exe)):
        toolsD['diamond'] = True
    if(which(hmmer_exe1) and which(hmmer_exe2)):
        toolsD['hmmer'] = True
    if(which(salmon_exe)):
        toolsD['salmon'] = True
    return toolsD


def check_dir_and_log():
    if(not os.path.isdir(PATH_TOOLS)):
        os.mkdir(PATH_TOOLS)
    sys.path.append(PATH_TOOLS) 
    if(not os.path.isfile(os.path.join(PATH_TOOLS, 'tools_log'))):
        write_log({})

def print_install_instructions(toolsDt):
    instructions = "some instructions will go here"
    # how best to print?
    print(instructions)

def run_tasks(tasks, cpu=4):
    for t in tasks:
        print(t.name)
        t.stdout = os.path.join(PATH_TOOLS, t.name+'.stdout')
        t.stderr = os.path.join(PATH_TOOLS, t.name+'.stderr')

    s = Supervisor(tasks=tasks, force_run=False, log=tool_supervisor_log, cpu=cpu)
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
    tfile.extractall(PATH_TOOLS)
    os.remove(target+'.tar.gz')


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
    log = open(os.path.join(PATH_TOOLS, 'tool_log'), 'w')
    json.dump(log_table, log, sort_keys=True, indent=4)
    log.close()


def main(install=False, toolList = [], tool_check=True, cpu=4):
    check_dir_and_log()
    toolsD = {}
    for tool in toolList:
        toolsD[tool.lower()] = False
    if tool_check:
        toolsD = check_tools(toolsD) #return altered tools dictionary
    log_table = read_log()
    tasks = []
    partial_get = lambda a, b, c : get(log_table, a, b ,c, tool_check)
    if(install and platform.system().lower() == 'linux'):
        if(not toolsD.get('trinity', True)): # if the tool is not in dictionary, we don't want to install, so default = True (pretend it's there)
	    partial_get('tar', trinity_linux_url,  trinity_linux_target)
            install_trinity = install_trinity_task( trinity_linux_target, trinity_exe, [], False)
	    tasks.append(install_trinity)
	if(not toolsD.get('prinseq', True)):
	    partial_get('tar', prinseq_url, prinseq_target)
	    install_prinseq = install_prinseq_task(prinseq_target, prinseq_exe, [], False)
            tasks.append(install_prinseq)
	if(not toolsD.get('transdecoder', True)):
	    partial_get('tar', transdecoder_url, transdecoder_target)
	    install_transdecoder = install_transdecoder_task(transdecoder_target, transdecoder_exe, [], False)
	    tasks.append(install_transdecoder)
	if(not toolsD.get('transrate', True)):
	    partial_get('tar', transrate_linux_url, transrate_linux_target)
	    install_transrate = install_transrate_task(transrate_linux_target, transrate_exe, [], False)
	    tasks.append(install_transrate)
	if(not toolsD.get('busco', True)):
	    partial_get('tar', busco_url, busco_target)
	    install_busco = install_busco_task(busco_target, busco_exe, [], False)
	    tasks.append(install_busco)
	if(not toolsD.get('diamond', True)):
	    partial_get('tar', diamond_linux_url, diamond_linux_target)
	    install_diamond = install_diamond_task(diamond_linux_target, diamond_exe, [], False)
	    tasks.append(install_diamond)
	if(not toolsD.get('hmmer', True)):
	    partial_get('tar', hmmer_linux_url, hmmer_linux_target)
	    install_hmmer = install_hmmer_task(hmmer_linux_target, hmmer_exe1, hmmer_exe2, [], False)
	    tasks.append(install_hmmer)
	if(not toolsD.get('salmon', True)):
	    partial_get('tar', salmon_linux_url, salmon_linux_target)
	    install_salmon = install_salmon_task(salmon_linux_target, salmon_exe, [], False)
	    tasks.append(install_salmon)
#        if(not toolsD['trimmomatic']):
#	    download_external_tool(prinseq_url, trimmomatic_url_target)
#           install_trimmomatic = install_trimmomatic_task(trimmomatic_url_target, trimmomatic_exe, [], False)
#	    tasks.append(install_trimmomatic)
    else:
        for tool, val in toolsD.items():
	    if not val:
	        print_install_instructions(toolsD)
    run_tasks(tasks, cpu)
    write_log(log_table)


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument('--install', action='store_true', default=False)
    parser.add_argument('--hard', action='store_true', default=False)
    parser.add_argument('-t', '--tool', action='append', default=[])
    parser.add_argument('--cpu', type=int, default=4)
    args = parser.parse_args()
    main(args.install, args.tool, not args.hard, args.cpu)

