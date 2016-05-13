# author: bluegenes

from time import strftime
import tarfile
import zipfile
from mmt_defaults import PATH_TOOLS, PATH_ROOT
from external_tools import TOOLS_DICT
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

tool_supervisor_log = '{0!s}/.tool_supervisor_log'.format(PATH_TOOLS)

def tool_check(t,fullpaths_exe, exe, allow_path=False):
    external_tools = False
    path_var = False
    if all(exists(x) for x in fullpaths_exe):
        external_tools = True
    if allow_path:
        if all(which(x) for x in exe):
            path_var = True
	elif all(which(os.path.basename(x)) for x in exe): # is this necessary?
	    path_var = True
	    t.exe = [os.path.basename(x) for x in exe]
    return any([external_tools, path_var]) # return both vals to distinguish btwn our installs vs path installs...

def check_tools(toolsD):
    toolsToInstall = {}
    for name, t in toolsD.items():
        if not tool_check(t,t.full_exe,t.exe):
            toolsToInstall[name] = t
    return toolsToInstall

def check_dir_and_log():
    if(not os.path.isdir(PATH_TOOLS)):
        os.mkdir(PATH_TOOLS)
    if(not os.path.isfile(join(PATH_TOOLS,'tools_log'))):
        write_log({})

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

def safe_retrieve(source, target, urltype):
    extension = ''
    if all( [len(urltype) >0, not urltype.startswith('.')]):
       extension = '.' + urltype
    temp = target+'.temp'
    print('getting ' + source)
    try:
        urlretrieve(source, temp)
        os.rename(temp, target+extension)
    except:
        print("WARNING : could not retrieve " + source)
    try:        
        expand_target(target, extension)
    except:
        print("WARNING : could not extract " + target+extension)

def expand_target(target, extension):
   if (extension == '.tar.gz' or extension == '.tgz'):
        tfile = tarfile.open((target+extension), 'r:gz')
        tfile.extractall(PATH_TOOLS)
        os.remove(target+extension)
   elif (extension == '.gz'):
        f = gzip.open((target+extension), 'rb')
        g = open(target, 'wb')
        for line in f:
            g.write(line)
        f.close()
        g.close()
        os.remove(target+extension)
   elif (extension == '.zip'):
        z = zipfile.ZipFile((target+extension))
        z.extractall(PATH_TOOLS)
        os.remove(target+extension)
   else:
        raise ValueError('Can\'t expand '+target+ ' -- please check the urltype/extension. Acceptable values are: "", "zip", "gz", "tgz", or  "tar.gz".')

def get(log_table, urltype, source, target, file_check=True):
    if(file_check and exists(target)):
        return
    safe_retrieve(source, target, urltype)
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
    for toolname in toolList:
        t = TOOLS_DICT[toolname]
        toolsD[t.name] = t
    if tool_check: # --hard option == install all tools that were passed in, no matter what we already have
        toolsD = check_tools(toolsD) #dict with only tools that need to be installed
    log_table = read_log()
    tasks = []
    partial_get = lambda a, b, c : get(log_table, a, b ,c, tool_check)
    if(install and platform.system().lower() == 'linux'):
        for name, tool in toolsD.items():
            if tool.install: # don't download tools that are not openly licensed or do not provide linux binaries
                partial_get(tool.urltype, tool.url, tool.target)
                install_task = tool.install_task
                if install_task is not None:
                    tasks.append(install_task)
                for flag, exe in zip(tool.executeable_flags, tool.full_exe):
                    if(flag):
                        cmd = 'chmod u+x {0!s}'.format(exe)
                        subprocess.call(cmd, shell=True)
            else:
                print(tool.instructions)
    else:
        for name, tool in toolsD.items():
            print('\n Installation instructions for: ' + tool.name)
            print('\n\t Download tool at this link: ' + tool.url)
            print(tool.instructions)
    run_tasks(tasks, cpu)
    write_log(log_table)


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument('--install', action='store_true', default=False)
    parser.add_argument('--hard', action='store_true', default=False)
    parser.add_argument('-t', '--tool', action='append', default=[])
    parser.add_argument('--cpu', type=int, default=4)
    args = parser.parse_args()
    #if args.hard:
    #    args.install = True
    main(args.install, args.tool, not args.hard, args.cpu)

