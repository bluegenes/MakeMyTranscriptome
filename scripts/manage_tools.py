from time import strftime
import tarfile
from task_functions_v2 import (PATH_ROOT, PATH_TOOLS,
install_trinity_task, install_trimmomatic_task, install_prinseq_task,
install_transdecoder_task, install_hmmer_task, install_salmon_task,
install_busco_task, install_transrate_task)
import os
from os.path import exists, join
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


class tool_variables:

    trinity_linux_url = 'https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.1.1.tar.gz'
    trinity_linux_target = join(PATH_TOOLS, 'trinityrnaseq-2.1.1')
    trinity_exe = 'Trinity'

    trimmomatic_url = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip'
    trimmomatic_trinity_target = join(which('Trinity'), 'trinity-plugins/Trimmomatic')
    trimmomatic_url_target = join(PATH_TOOLS, 'Trimmomatic-0.35')
    trimmomatic_exe = 'Trimmomatic-0.35.jar'

    prinseq_url = 'http://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz'
    prinseq_target = join(PATH_TOOLS, 'prinseq-lite-0.20.4')
    prinseq_exe = 'prinseq-lite.pl'

    transdecoder_url = 'https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz'
    transdecoder_target = join(PATH_TOOLS, 'TransDecoder-2.0.1')
    transdecoder_exe1 = 'TransDecoder.Predict'
    transdecoder_exe2 = 'TransDecoder.LongOrfs'

    transrate_linux_url = 'https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz'
    transrate_linux_target = join(PATH_TOOLS, 'transrate-1.0.1-linux-x86_64')
    transrate_exe = 'transrate'

    busco_url = 'http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz'
    busco_target = join(PATH_TOOLS, 'BUSCO_v1.1b1')
    busco_exe = 'BUSCO_v1.1b1.py'

    hmmer_linux_url = 'http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz'
    hmmer_linux_target = join(PATH_TOOLS, 'hmmer-3.1b2-linux-intel-x86_64')
    hmmer_exe1 = 'hmmscan' 
    hmmer_exe2 = 'hmmpress'

    diamond_linux_url = 'http://github.com/bbuchfink/diamond/releases/download/v0.7.10/diamond-linux64.tar.gz'
    diamond_linux_target = join(PATH_TOOLS, '_diamond')
    diamond_exe = 'diamond'

    #salmon_linux_url = 'https://github.com/COMBINE-lab/salmon/releases/download/v0.5.1/SalmonBeta-0.5.1_DebianSqueeze.tar.gz'
    #salmon_linux_target = join(PATH_TOOLS, 'SalmonBeta-0.5.1_DebianSqueeze')
    salmon_linux_url = 'https://github.com/COMBINE-lab/salmon/releases/download/v0.6.0/SalmonBeta-0.6.0_DebianSqueeze.tar.gz'
    salmon_linux_target = join(PATH_TOOLS, 'SalmonBeta-0.6.1_DebianSqueeze')
    salmon_exe = 'salmon'

    # trinity_version_check = 'Trinity --version'
    # prinseq_version_check = 'perl prinseq-lite.pl -version'
    # transrateversion_check = 'transrate --version'

tv = tool_variables
tool_supervisor_log = '{0!s}/.tool_supervisor_log'.format(PATH_TOOLS)


def tools_join(p):
    return join(PATH_TOOLS, p)


def check_tools(toolsD):
    if(which(tv.trinity_exe) or exists(tools_join(tv.trinity_exe))):
        toolsD['trinity'] = True
    if(which(tv.prinseq_exe) or exists(tools_join(tv.prinseq_exe))):
        toolsD['prinseq'] = True
    if(which(tv.transrate_exe or exists(tools_join(tv.transrate_exe)))):
        toolsD['transrate'] = True
    if(which(tv.busco_exe) or exists(tools_join(tv.busco_exe))):
        toolsD['busco'] = True
    if(which(tv.salmon_exe) or exists(tools_join(tv.salmon_exe))):
        toolsD['salmon'] = True
    if(which(tv.diamond_exe) or exists(tools_join(tv.diamond_exe))):
        toolsD['diamond'] = True
    if(all([which(tv.transdecoder_exe1), which(tv.transdecoder_exe2)]) or
       all([exists(tools_join(tv.transdecoder_exe1)), exists(tools_join(tv.transdecoder_exe2))])):
        toolsD['transdecoder'] = True
    if(all([which(tv.transdecoder_exe1), which(tv.transdecoder_exe2)]) or
       all([exists(tools_join(tv.transdecoder_exe1)), exists(tools_join(tv.transdecoder_exe2))])):
        toolsD['transdecoder'] = True
    if(all([which(tv.hmmer_exe1), which(tv.hmmer_exe2)]) or
       all([exists(tools_join(tv.hmmer_exe1)), exists(tools_join(hmmer_exe2))])):
        toolsD['hmmer'] = True
    return toolsD


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
        toolsD[tool.lower()] = False
    if tool_check:
        toolsD = check_tools(toolsD) #return altered tools dictionary
    else:
        install = True
    log_table = read_log()
    tasks = []
    partial_get = lambda a, b, c : get(log_table, a, b ,c, tool_check)
    if(install and platform.system().lower() == 'linux'):
        if('trinity' in toolsD):  # if the tool is not in dictionary, we don't want to install
            install_trinity = install_trinity_task( tv.trinity_linux_target, tv.trinity_exe, [], False)
            partial_get('tar', tv.trinity_linux_url,  tv.trinity_linux_target)
            tasks.append(install_trinity)
        if('prinseq' in toolsD):
            partial_get('tar', tv.prinseq_url, tv.prinseq_target)
            install_prinseq = install_prinseq_task(tv.prinseq_target, tv.prinseq_exe, [], False)
            tasks.append(install_prinseq)
        if('transdecoder' in toolsD):
            partial_get('tar', tv.transdecoder_url, tv.transdecoder_target)
            install_transdecoder = install_transdecoder_task(tv.transdecoder_target, tv.transdecoder_exe1, tv.transdecoder_exe2, [], False)
            tasks.append(install_transdecoder)
        if('transrate' in toolsD):
            partial_get('tar', tv.transrate_linux_url, tv.transrate_linux_target)
            install_transrate = install_transrate_task(tv.transrate_linux_target, tv.transrate_exe, [], False)
            tasks.append(install_transrate)
        if('busco' in toolsD):
            partial_get('tar', tv.busco_url, tv.busco_target)
            install_busco = install_busco_task(tv.busco_target, tv.busco_exe, [], False)
            tasks.append(install_busco)
        if('diamond' in toolsD):
            partial_get('tar', tv.diamond_linux_url, tv.diamond_linux_target)
        if('hmmer' in toolsD):
            partial_get('tar', tv.hmmer_linux_url, tv.hmmer_linux_target)
            install_hmmer = install_hmmer_task(tv.hmmer_linux_target, tv.hmmer_exe1, tv.hmmer_exe2, [], False)
            tasks.append(install_hmmer)
        if('salmon' in toolsD):
            partial_get('tar', tv.salmon_linux_url, tv.salmon_linux_target)
            install_salmon = install_salmon_task(tv.salmon_linux_target, tv.salmon_exe, [], False)
            tasks.append(install_salmon)
#        if(not toolsD['trimmomatic']):
#        download_external_tool(prinseq_url, trimmomatic_url_target)
#           install_trimmomatic = install_trimmomatic_task(trimmomatic_url_target, trimmomatic_exe, [], False)
#        tasks.append(install_trimmomatic)
    else:
        for tool, val in toolsD.items():
            if(not val):
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

