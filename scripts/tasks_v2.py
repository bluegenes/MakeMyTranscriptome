import os
import sys
import time
import signal
import subprocess
import pickle
from itertools import chain
import warnings


def time_to_hms(delta):
    m, s = divmod(delta, 60)
    h, m = divmod(m, 60)
    return (h, m, s)


class Task:

    class ExitCodeException(Exception):
        exit_code = None

    class TaskException(Exception):
        pass

    def __init__(self, command, dependencies=[], targets=[], cpu=1, name='Anonymous_Task', stderr=None, stdout=None, error_check=None, max_wall_time=float('inf')):
        try:
            if(stderr is not None):
                f = open(stderr, 'a')
                f.close()
            if(stdout is not None):
                f = open(stdout, 'a')
                f.close()
        except:
            raise
        self.command = command
        self.dependencies = dependencies
        self.cpu = cpu
        self.name = name
        self.stdout = stdout
        self.stderr = stderr
        self.targets = targets
        self.error_check = error_check if(error_check!=None) else lambda exit_code : exit_code!=0
        self.opened_files = []
        self.process = None
        self.exit_code = None
        self.soft_finished_status = False
        self.start_time = None

    def checkDependencies(self):
        for d in self.dependencies:
            if( isinstance(d, Task) or isinstance(d, Supervisor)):
                try:
                    if( not d.finished()):
                        return False
                except (Task.ExitCodeException, Task.TaskException):
                    return False
            elif( callable(d) ): 
                if( not d() ):
                    return False
            else:
                error_message = 'Unable to check depencies of task '+self.name
                error_message+='. \nTask dependencies must not conatin anything that is not a function or a Task. '
                raise self.TaskException(error_message)
        return True

    def start(self):
        if(self.stderr!=None):
            err = open(self.stderr,'w',1)
            err.write(str(self.command)+'\n\n')
            self.opened_files.append(err)
        else:
            err = None
        if(self.stdout!=None):
            out = open(self.stdout,'w',1)
            out.write(str(self.command)+'\n\n')
            self.opened_files.append(out)
        else:
            out = None
        self.start = time.time()
        temp = subprocess.Popen(self.command,shell=True,stdout=out,stderr=err)
        self.process = temp

    def finished(self):
        if(self.soft_finished_status):
            return True
        if(self.start_time is not None and float(time.time()-self.start_time)/60 > self.max_wall_time):
            err_mess = ('Task {0!s} has been running for greater than its maximum wall time, '
                        '{1!s}m, and has been aborted. This is likely an external error and '
                        'trying again is recommended.').format(self.name, self.max_wall_time)
            self.killRun()
            raise self.TaskException(err_mess)
        if(self.process==None):
            return False
        self.process.poll()
        exit_code = self.process.returncode
        if(exit_code == None):
            return False
        else:
            self.exit_code=exit_code
            self.opened_files = [f.close() for f in self.opened_files]
            self.opened_files = []
            if(self.error_check(exit_code)):
                error_message = 'Task '+ self.name+' seems to have failed with exit_code '+str(exit_code)+'.'
                if(self.stderr!=None):
                    error_message+=' Consult '+self.stderr+' for more info.' 
                inst = self.ExitCodeException(error_message)
                inst.exit_code = exit_code
                raise inst
            for t in self.targets:
                if(not os.path.exists(t)):
                    error_message = 'Task '+self.name+' encountered an unexpected error resulting in it failing'
                    error_message+= ' to produce its target files.'
                    if(self.stderr!=None):
                        error_message+=' Consult '+self.stderr+' for more info.'
                    raise self.TaskException(error_message)
            return True

    def run(self):
        self.start()
        self.process.wait()
        self.finished()

    def killRun(self):
        try:
            self.process.kill()
        except:
            pass
        for f in self.opened_files:
            f.close()
        self.opened_files = []

    def skipable(self, history):
        if(self.soft_finished_status):
            return True
        if(self.command not in history):
            return False
        for t in self.dependencies:
            if(isinstance(t, Task) or isinstance(t, Supervisor)):
                if(not t.skipable(history)):
                    return False
        for t in self.targets:
            if(not os.path.exists(t)):
                return False
        self.soft_finished_status = True
        return True

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.__str__()


class Supervisor:

    STATE_INITIALIZED = 'initialized'
    STATE_ERR = 'failed'
    STATE_SKIPPED = 'skipped'
    STATE_FINISHED = 'executed'
    STATE_RUNNING = 'started'
    STATE_REMOVED = 'removed'
    history_log_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '.task_log')


    def __init__(self, tasks=[], dependencies=[], cpu=float('inf'), name='Supervisor', delay=1, force_run=False, email=None, email_interval=30, log=None):
        self.cpu = cpu
        self.name = name
        self.delay = delay
        self.dependencies = dependencies
        self.force_run = force_run
        self.email = email
        self.email_interval = email_interval
        self.last_email = time.time()
        self.log_path = log if(log is not None) else name+'.run_log'
        self.log_str = ''
        self.task_map = {}
        self.task_status = {}
        self.errors = []
        self.targets = []
        self.tasks = set()
        for t in tasks:
            self.add_task(t)

    def run(self):
        self.log_file = open(self.log_path,'w',1)
        self.last_email = time.time()
        self.tasks_to_run = set([t for t in self.tasks])
        self.tasks_running = set()
        cur_cpu = 0
        history = self.get_history()
        history_update = set()
        signal.signal(signal.SIGTERM, lambda *x: self.killRun())
        try:
            # execute run
            while(len(self.tasks_to_run) > 0 or len(self.tasks_running) > 0):
                size_tasks_to_run = len(self.tasks_to_run)
                # handling finished tasks
                for t in [task for task in self.tasks_running]:
                    temp = self.task_status[t]
                    try:
                        if(t.finished()):
                            self.tasks_running.remove(t)
                            cur_cpu -= t.cpu
                            temp['stop'] = int(time.time())
                            temp['state'] = Supervisor.STATE_FINISHED
                            temp['exit_code'] = t.exit_code
                            h, m, s = time_to_hms(temp['stop'] - temp['start'])
                            temp['message'] = 'Completed in {0!s}h {1!s}m {2!s}s'.format(h, m, s)
                            history_update.add(t.command)
                            self.log(t.name+':'+self.task_status[t]['state']+':'+time.asctime()+'\n\n')
                    except (Task.ExitCodeException, Task.TaskException) as inst:
                        cur_cpu -= t.cpu
                        self.tasks_running.remove(t)
                        self.errors.append(inst)
                        temp['stop'] = int(time.time())
                        temp['state'] = Supervisor.STATE_ERR
                        temp['exit_code'] = t.exit_code
                        h, m, s = time_to_hms(temp['stop'] - temp['start'])
                        temp['message'] = 'Failed in {0!s}h {1!s}m {2!s}s'.format(h, m, s)
                        self.log(t.name+':'+self.task_status[t]['state']+':'+time.asctime()+'\n\n')
                        self.__removeTaskPath__(t)
                # starting execution of tasks
                for t in [task for task in self.tasks_to_run]:
                    temp = self.task_status[t]
                    if(t.cpu+cur_cpu > self.cpu):
                        continue
                    elif(not t.checkDependencies()):
                        continue
                    elif(not self.force_run and t.skipable(history)):
                        temp['state'] = Supervisor.STATE_SKIPPED
                        self.log(t.name+':'+temp['state']+'\n')
                        self.tasks_to_run.remove(t)
                        continue
                    else:
                        cur_cpu += t.cpu
                        self.tasks_running.add(t)
                        self.tasks_to_run.remove(t)
                        t.start()
                        temp['state'] = Supervisor.STATE_RUNNING
                        temp['start'] = int(time.time())
                        self.log(t.name+':'+self.task_status[t]['state']+':'+time.asctime()+'\n\n')
                        self.log_file.write(t.command+'\n\n')
                tasks_to_run_delta = size_tasks_to_run - len(self.tasks_to_run)
                if(len(self.tasks_running) == 0 and tasks_to_run_delta == 0):
                    break
                time.sleep(self.delay)
            # handle all errors
            if(self.errors != [] or (len(self.tasks_running) == 0 and len(self.tasks_to_run) != 0)):
                err_str = '\n\n'
                if(len(self.tasks_running) == 0 and len(self.tasks_to_run) != 0):
                    err_str += 'Unable to resolve dependencies during execution of '+self.name
                    err_str += '. The following tasks could not be executed:\n'
                    err_str += '\n'.join(['\t'+t.name for t in self.tasks_to_run])
                if(self.errors != []):
                    err_str += '\nEncountered an unexpected Error in the following tasks:\n'
                    for t in self.task_status:
                        if(self.task_status[t]['state'] == Supervisor.STATE_ERR):
                            err_str+= '\tName - {0!s} : Message - {1!s} : Exit Code - {2!s}\n'.format(
                                t.name,self.task_status[t]['message'],self.task_status[t]['exit_code'])
                    removed_tasks = []
                    for t in self.task_status:
                        if(self.task_status[t]['state']==Supervisor.STATE_REMOVED):
                            removed_tasks.append('\tName - {0!s} : Message - {1!s} : Exit Code- {2!s}\n'.format(
                                t.name,self.task_status[t]['message'],self.task_status[t]['exit_code']))
                    if(len(removed_tasks)>0):
                        err_str+='\nAs a result of the above errors, the following tasks could not be executed:\n'
                        err_str+=''.join(removed_tasks)
                    err_str+='\nErrors Reported:\n'
                    for e in self.errors:
                        err_str+=str(e)+'\n'
                raise Exception(err_str)            
        except BaseException as inst:
            self.killRun()
            self.log_file.write(str(inst))
            raise
        finally:
            self.write_history(history_update)
            self.send_email('', subject='Pipeline Finished')

    def __removeTaskPath__(self, task):
        flag = True
        removed = set([task])
        check_intersection = lambda s, l : any(e in s for e in l)
        while(flag):
            flag = False
            for t in [t for t in self.tasks_to_run]:
                sup_deps = [d for d in t.dependencies if(isinstance(d, Supervisor))]
                sup_checks = [check_intersection(removed, s.tasks) for s in sup_deps]
                if(check_intersection(removed, t.dependencies) or any(sup_checks)):
                    temp = self.task_status[t]
                    self.tasks_to_run.remove(t)
                    removed.add(t)
                    temp['message'] = 'Never Started'
                    temp['state'] = Supervisor.STATE_REMOVED
                    flag=True

    def get_history(self):
        if(os.path.isfile(self.history_log_path)):
            history_file = open(self.history_log_path, 'rb')
            ret = pickle.load(history_file)
            history_file.close()
        else:
            ret = set()
        return ret

    def write_history(self, update):
        history = self.get_history()
        history = history.union(update)
        history_file = open(self.history_log_path, 'wb')
        pickle.dump(history, history_file)
        history_file.close()

    def log(self, message):
        self.log_str += message
        self.log_file.write(message)
        print(message)
        if(time.time() > self.last_email + self.email_interval):
            self.last_email = time.time()
            self.send_email(self.log_str, 'Pipeline Running Update')
            self.log_str = ''

    def finished(self):
        for t in self.tasks:
            if(not t.finished()):
                return False
        return True

    def skipable(self, history):
        for t in self.tasks:
            if(not t.skipable(history)):
                return False
        return True

    def send_email(self, message, subject=''):
        if(self.email is None):
            return
        else:
            # message = ''.join([c if(c!='\n') else '\\n'for c in message])
            cmd = "echo '{0!s}' | mail -s '{1!s}' '{2!s}'".format(message, subject, self.email)
            subprocess.call(cmd, shell=True)

    def killRun(self):
        for t in self.tasks_running:
            t.killRun()
        self.running = []

    def add_task(self, task):
        if(isinstance(task, Supervisor)):
            for t in task.tasks:
                t.dependencies.extend(task.dependencies)
                self.add_task(t)
        if(isinstance(task, Task)):
            if(task.cpu > self.cpu):
                err_mess = ('Task {0!s} has a higher cpu than this supervisor, {1!s}.\nYou '
                            'must increase this supervisor\'s cpu or decrease {0!s}\'s cpu.')
                err_mess = err_mess.format(task.name, self.name)
                raise Exception(err_mess)
            if(task.name in self.task_map):
                warnings.warn('A task named '+task.name+' already exists.')
            self.task_map[task.name] = task
            self.tasks.add(task)
            self.task_status[task] = {'state': Supervisor.STATE_INITIALIZED, 'exit_code': None,
                                      'message': None, 'start': None, 'stop': None}
            self.targets.extend(task.targets)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.__str__()



if(__name__=='__main__'):
    t1 = Task('py3 ..\\..\\test.py 4.4',dependencies=[],name='t1')
    t2 = Task('py3 ..\\..\\test.py 6',dependencies=[t1],name='t2')
    def t2_dep():
        t2.command = 'py3 ..\\..\\test.py 8'
        return True
    t2.dependencies.append(t2_dep)
    t3 = Task('py3 ..\\..\\test.py 3',dependencies=[],name='t3')
    t4 = Task('py3 ..\\..\\test.py 2',dependencies=[t3],name='t4')
    t5 = Task('py3 ..\\..\\test.py 1',dependencies=[],name='t5')
    t6 = Task('py3 ..\\..\\test.py 7',dependencies=[],name='t6')
    s1 = Supervisor([t1,t2],name='s1')
    s2 = Supervisor([t3,t4],name='s2',dependencies=[s1])
    s3 = Supervisor([s1,s2],name='s3')
    s4 = Supervisor([t5,t6],name='s4',dependencies=[t1])
    s = Supervisor([s3,s4],force_run=False)
    s.run()


