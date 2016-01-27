# author: bluegenes

import os
from os.path import join, exists
from tasks_v2 import Task

class ext_tool:
    def __init__(self,name,url,target,executables,instructions,urltype='tar.gz',dependencies = []):
        self.name = name
        self.url = url
        self.target = target
        self.exe = executables
        self.target = target
        self.instructions = instructions
        self.urltype = urltype
	self.dependencies = dependencies
        self.full_exe = [join(self.target, x) for x in self.exe]
        self.install_task = None
        self.check_install_task = None
	
    def set_install(self, cmd, task_deps=[], log_flag=True):
        install_trgs = self.full_exe
        cd_cmd = 'cd {0!s}; '.format(self.target)
        install_cmd = cd_cmd + cmd
        install_name ='install_' + self.name
	out,err = (None,None)
#        out, err = GEN_LOGS(name) if(log_flag) else (None, None)
        self.install_task=Task(command=install_cmd,dependencies=task_deps,targets=install_trgs,name=install_name,stdout=out,stderr=err)

    def check_install(self,check_cmd,task_deps=[], log_flag=True):
        install_trgs = self.full_exe
        check_name ='check_install_' + self.name
	out,err = (None,None)
 #       out, err = GEN_LOGS(name) if(log_flag) else (None, None)
        self.check_install_task=Task(command=check_cmd,dependencies=task_deps,targets=install_trgs,name=check_install_name,stdout=out,stderr=err) 

    def change_exe_fullpath(self, path):
        self.full_exe = [join(path, x) for x in self.exe]
    
    def __call__(self):
        return [self.name, self.url, self.target, self.exe,self.urltype, self.install, self.instructions]


