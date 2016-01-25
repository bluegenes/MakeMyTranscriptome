from time import strftime
import tarfile
from task_functions_v2 import (PATH_ROOT, PATH_TOOLS,
install_trinity_task, install_trimmomatic_task, install_prinseq_task,
install_transdecoder_task, install_hmmer_task, install_salmon_task,
install_busco_task, install_transrate_task)
import os
from os.path import exists, join, split
import json
import argparse
import sys
import functools
import subprocess
import platform
from tasks_v2 import Supervisor
if(sys.version[0] == '3'):
    from urllib.request import urlretrieve, ContentTooShortError
    from shutil import which as which
else:
    from urllib import urlretrieve, ContentTooShortError
    from  py2_which import which_python2 as which

class tool_variables:

    def __init__(self, name, url, target, executables, install_task):
        self.name = name
        self.url = url
        self.targets = target
        self.exe = executables
        self.target = target
	if install_task is not None:
	    self.install = install_task(self.target, self.exe, [], False)
	else:
	    self.install = None

    def add_exe(self, exe):
        self.executables.append(exe)

    def __call__(self):
        return [self.name, self.url, self.target, self.executables, self.install]

########### links and targets ###########
trinity_linux_url = 'https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz'
trinity_linux_target = join(PATH_TOOLS, 'trinityrnaseq-2.1.1')
trinity_exe = ['Trinity']

trimmomatic_url = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip'
trimmomatic_target = join(PATH_TOOLS, 'Trimmomatic-0.35')
trimmomatic_exe = ['Trimmomatic-0.35.jar']
#trimmomatic_trinity_target = join(split(which('Trinity'))[0], 'trinity-plugins/Trimmomatic')

prinseq_url = 'http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz'
prinseq_target = join(PATH_TOOLS, 'prinseq-lite-0.20.4')
prinseq_exe = ['prinseq-lite.pl']

transdecoder_url = 'https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz'
transdecoder_target = join(PATH_TOOLS, 'TransDecoder-2.0.1')
transdecoder_exe = ['TransDecoder.Predict','TransDecoder.LongOrfs']

transrate_linux_url = 'https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz'
transrate_linux_target = join(PATH_TOOLS, 'transrate-1.0.1-linux-x86_64')
transrate_exe = ['transrate']

busco_url = 'http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz'
busco_target = join(PATH_TOOLS, 'BUSCO_v1.1b1')
busco_exe = ['BUSCO_v1.1b1.py']

hmmer_linux_url = 'http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz'
hmmer_linux_target = join(PATH_TOOLS, 'hmmer-3.1b2-linux-intel-x86_64')
hmmer_exe = ['hmmscan','hmmpress']

diamond_linux_url = 'http://github.com/bbuchfink/diamond/releases/download/v0.7.10/diamond-linux64.tar.gz'
diamond_target = join(PATH_TOOLS, '_diamond')
diamond_exe = ['diamond']

salmon_linux_url = 'https://github.com/COMBINE-lab/salmon/releases/download/v0.6.0/SalmonBeta-0.6.0_DebianSqueeze.tar.gz'
salmon_linux_target = join(PATH_TOOLS, 'SalmonBeta-0.6.1_DebianSqueeze')
salmon_exe = ['salmon']

############# define the tools ##################
tv = tool_variables
trinity_tool = tv('trinity', trinity_linux_url, trinity_linux_target, trinity_exe, install_trinity_task)
trimmomatic_tool = tv('trimmomatic', trimmomatic_url, trimmomatic_target, trimmomatic_exe, install_trimmomatic_task)
prinseq_tool = tv('prinseq', prinseq_url, prinseq_target, prinseq_exe, install_prinseq_task) 
transdecoder_tool = tv('transdecoder', transdecoder_url, transdecoder_target, transdecoder_exe, install_transdecoder_task)
transrate_tool = tv('transrate', transrate_linux_url, transrate_linux_target, transrate_exe, install_transrate_task)
busco_tool = tv('busco', busco_url, busco_target, busco_exe, install_busco_task)
hmmer_tool = tv('hmmer',hmmer_linux_url, hmmer_linux_target, hmmer_exe, install_hmmer_task)
diamond_tool = tv('diamond', diamond_linux_url, diamond_target, diamond_exe, None)
salmon_tool = tv('salmon', salmon_linux_url, salmon_linux_target, salmon_exe, install_salmon_task)

tool_list = [trinity_tool, trimmomatic_tool, prinseq_tool, transdecoder_tool, transrate_tool,busco_tool, hmmer_tool, diamond_tool, salmon_tool]
tvTools = {tool.name:tool for tool in tool_list}
##################################################

tool_supervisor_log = '{0!s}/.tool_supervisor_log'.format(PATH_TOOLS)

def tools_join(p):
    return join(PATH_TOOLS, p)

def tool_check(exe_list, allow_path=False):
    external_tools = False
    path_var = False
    et_paths = [tools_join(x) for x in exe_list]
    if all(exists(x) for x in et_paths):
	external_tools = True
    if allow_path:
        if all(which(x) for x in exe_list):
	    path_var = True
    return any( [external_tools, path_var]) # just remove the 'any' and return both vals if want to be able to distinguish btwn our installs vs theirs...

def check_tools(toolsD):
    toolsToInstall = {}
    for name, t in toolsD.items():
        if not tool_check(t.exe):
            toolsToInstall[name] = t
    return toolsToInstall

def check_dir_and_log():
    if(not os.path.isdir(PATH_TOOLS)):
        os.mkdir(PATH_TOOLS)
    if(not os.path.isfile(tools_join('tools_log'))):
        write_log({})

def print_install_instructions(toolsDt):
    instructions = "some instructions will go here"
    # how best to print?
    print(instructions)


def run_tasks(tasks, cpu=4):
    for t in tasks:
        print(t.name)
        t.stdout = join(PATH_TOOLS, t.name+'.stdout')
        t.stderr = join(PATH_TOOLS, t.name+'.stderr')

    s = Supervisor(tasks=tasks, force_run=False, log=tool_supervisor_log, cpu=cpu)
    s.run()
    for t in tasks:  #if everything executes properly, rm the task logs
        if exists(t.stdout):
            os.remove(t.stdout)
        if exists(t.stderr):
            os.remove(t.stderr)


def safe_retrieve(source, target):
    temp = target+'.temp'
    print('getting ' + source)
    urlretrieve(source, temp)
    os.rename(temp, target)


def url_unzip(source, target):
    print('getting '+source)
    gztemp = target+'.gz'
    urlretrieve(source, gztemp)
    f = gzip.open(gztemp, 'rb')
    g = open(target, 'wb')
    for line in f:
        g.write(line)
    f.close()
    g.close()
    os.remove(gztemp)


def tar_retrieve(source, target):
    print('getting '+source)
    tartemp = target+'.tar.gz'
    urlretrieve(source, tartemp)
    # safe_retrieve(source, target+'.tar.gz')
    tfile = tarfile.open(tartemp, 'r:gz')
    tfile.extractall(PATH_TOOLS) #hmmm
    os.remove(tartemp)


def get(log_table, flag, source, target, file_check=True):
    if(file_check and exists(target)):
        return
    try:
        if(flag == 'gz'):
            url_unzip(source, target)
        elif(flag == ''):
            safe_retrieve(source, target)
        elif(flag == 'tar'):
            tar_retrieve(source, target)
        else:
            raise ValueError('The flag used must be "", "gz", or "tar".')
    except ContentTooShortError:
        print('failed to install {0!s}'.format(source))
    basename = os.path.basename(target)
    log_table[basename] = strftime('%b-%d-%Y')


def read_log():
    log = open(join(PATH_TOOLS, 'tool_log'))
    log_table = json.load(log)
    log.close()
    return log_table


def write_log(log_table):
    log = open(join(PATH_TOOLS, 'tool_log'), 'w')
    json.dump(log_table, log, sort_keys=True, indent=4)
    log.close()


def main(install=False, toolList = [], tool_check=True, cpu=4):
    check_dir_and_log()
    toolsD = {}
    for tool in toolList:
        t = tvTools[tool]
        toolsD[t.name] = t
	print(toolsD)
	print(t.name, t.url, t.exe, t.target, t.install)
    if tool_check: # --hard option --> don't want to check if installed
        toolsD = check_tools(toolsD) #return altered tools dictionary
    else:
        install = True
    log_table = read_log()
    tasks = []
    partial_get = lambda a, b, c : get(log_table, a, b ,c, tool_check)
    if(install and platform.system().lower() == 'linux'):
	for tool in toolD:
	    #name = tool.name
	    partial_get('tar', tool.url, tool.target)
	    install_task = tool.install(tool.target, tool.exe, [], False)
	    tasks.append(install_task)

#if('trinity' in toolsD):  # if the tool is not in dictionary, we don't want to install
#            install_trinity = install_trinity_task( tv.trinity_linux_target, tv.trinity_exe, [], False)
#            partial_get('tar', tv.trinity_linux_url,  tv.trinity_linux_target)
#            tasks.append(install_trinity)
 #       if('prinseq' in toolsD):
 #           partial_get('tar', tv.prinseq_url, tv.prinseq_target)
 #           install_prinseq = install_prinseq_task(tv.prinseq_target, tv.prinseq_exe, [], False)
 #           tasks.append(install_prinseq)
 #       if('transdecoder' in toolsD):
 #           partial_get('tar', tv.transdecoder_url, tv.transdecoder_target)
#            install_transdecoder = install_transdecoder_task(tv.transdecoder_target, tv.transdecoder_exe1, tv.transdecoder_exe2, [], False)
#            tasks.append(install_transdecoder)
#        if('transrate' in toolsD):
#            partial_get('tar', tv.transrate_linux_url, tv.transrate_linux_target)
#            install_transrate = install_transrate_task(tv.transrate_linux_target, tv.transrate_exe, [], False)
#            tasks.append(install_transrate)
#        if('busco' in toolsD):
#            partial_get('tar', tv.busco_url, tv.busco_target)
#            install_busco = install_busco_task(tv.busco_target, tv.busco_exe, [], False)
#            tasks.append(install_busco)
#        if('diamond' in toolsD):
#            partial_get('tar', tv.diamond_linux_url, tv.diamond_linux_target)
#        if('hmmer' in toolsD):
#            partial_get('tar', tv.hmmer_linux_url, tv.hmmer_linux_target)
#            install_hmmer = install_hmmer_task(tv.hmmer_linux_target, tv.hmmer_exe1, tv.hmmer_exe2, [], False)
#            tasks.append(install_hmmer)
#        if('salmon' in toolsD):
#            partial_get('tar', tv.salmon_linux_url, tv.salmon_linux_target)
#            install_salmon = install_salmon_task(tv.salmon_linux_target, tv.salmon_exe, [], False)
#            tasks.append(install_salmon)
##        if('trimmomatic' in toolsD):
#            partial_get('tar', tv.trimmomatic_url, tv.trimmomatic_target)
#            install_trimmomatic = install_trimmomatic_task(tv.trimmomatic_target, tv.trimmomatic_exe, [], False)
#            tasks.append(install_trimmomatic)
    else:
        for tool in toolsD:
            print_install_instructions(tool)
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

