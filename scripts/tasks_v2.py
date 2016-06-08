import os
import sys
import time
import signal
import subprocess
import warnings
import platform
import shlex


def time_to_hms(delta):
    ''' Convert some time in seconds to a tuple of (hours, minutes, seconds).
    '''
    m, s = divmod(delta, 60)
    h, m = divmod(m, 60)
    return (h, m, s)


def check_foreground():
    ''' function that returns True if the process executing this module
        is in foreground. If in background, this module won't
        print to stderr/stdout.
    '''
    if(not platform.system().lower().startswith('linux')):
        return True
    try:
        if(os.getpgrp() == os.tcgetpgrp(sys.stdout.fileno())):
            return True
    except OSError:
        pass
    return False


def not_zero(i):
    ''' function that returns True if i is not zero. the funciton is equivalent
        to lambda i: i!=0
    '''
    return i != 0


class Task:
    ''' Object that acts as the basic unit of a pipeline. Command line jobs
        can be defined and exucted using a Task object
    '''

    class ExitCodeException(Exception):
        ''' An exception raised by a call to self.finished() if the subprocess
            exits with a bad exit code.
        '''
        exit_code = None

    class TaskException(Exception):
        ''' A more general TaskException that can be thrown be a task object
            and caught by outside funcitons.
        '''
        pass

    def __init__(
      self, command, dependencies=[], targets=[], cpu=1, name='unnamed_task',
      error_check=not_zero, max_wall_time=float("inf"), **popen_args):
        ''' The __init__ for the Task object has nine parameters described
            below. It also supports all keyword arguments of the
            subprocess.Popen object:
            args -
                almost equivalent to the subprocess.Popen argument args. If
                a string is given, it will be split using shlex.split before
                being passed to Popen
            dependencies -
                a list containing Task objects, Supervisor objects, or unit
                functions. A call to self.checkDependencies will return True
                if and only if all Task objects in dependencies are finished,
                all Supervisor objects in dependencies are finished, and all
                unit funcitons return True when called.
            targets -
                a list of strings that are the paths to the output files
                generated by running the task.
            cpu -
                the number of threads or cpu intensivity of the command. This
                value is used by Supervisor objects. Default=1
            name -
                A string that will be used for representing the task for more
                convienient logging. Default="unnamed_task"
            error_check -
                a function f(int)->Bool that the exit code
                from running command. The function should return True if the
                exit code implies an error was encountered while running
                command, else False. Default = lambda i: i != 0
            max_wall_time -
                the amount of time, in minutes, that the command should be
                allowed to run before being stopped by a call to
                self.finished(). Default = float('inf')
        '''
        self.command = command
        self.targets = targets
        self.dependencies = dependencies
        self.cpu = cpu
        self.name = name
        self.error_check = error_check
        self.max_wall_time = max_wall_time
        self.popen_args = popen_args
        self.process = None
        self.opened_files = []
        self.exit_code = None
        self.soft_finished_status = False
        self.start_time = -1

    def _handleFilePopenArgs(self, arg, type_flag='r'):
        ''' Allows for strings to be passed in popen args where files would
            normally be required. Its a helper function.
        '''
        if(arg in self.popen_args and isinstance(self.popen_args[arg], str)):
            temp_file = open(self.popen_args[arg], type_flag)
            self.opened_files.append(temp_file)
        elif(arg in self.popen_args):
            temp_file = self.popen_args[arg]
        else:
            temp_file = None
        return temp_file

    def start(self):
        ''' Method used to start the execution of this Task. self.command will
            be executed using the subprocess.Popen constructor with shell=True.
        '''
        print(self.name)
        stdin = self._handleFilePopenArgs("stdin")
        stdout = self._handleFilePopenArgs("stdout", "w")
        stderr = self._handleFilePopenArgs("stderr", "w")
        file_args = ["stdin", "stdout", "stderr"]
        args = {k: self.popen_args[k] for k in self.popen_args if(k not in file_args)}
        self.start_time = time.time()
        cmd = shlex.split(self.command) if(isinstance(self.command, str)) else self.command
        self.process = subprocess.Popen(cmd, stdin=stdin, stdout=stdout, stderr=stderr, **args)

    def checkDependencies(self):
        ''' Method will check all dependencies of this object.
            self.dependencies is a list containing Task objects, Supervisor
            objects, or unit functions. A call to self.checkDependencies will
            return True if and only if all Task objects in dependencies are
            finished, all Supervisor objects in dependencies are finished, and
            all unit funcitons return True when called.
        '''
        for d in self.dependencies:
            if(isinstance(d, Task) or isinstance(d, Supervisor)):
                try:
                    if(not d.finished()):
                        return False
                except (Task.ExitCodeException, Task.TaskException):
                    return False
            elif(callable(d)):
                if(not d()):
                    return False
            else:
                err_msg = ('Unable to check depencies of task {0!1}. \nTask '
                           'dependencies must not conatin anything that is not '
                           'a function or a Task.').format(self.name)
                raise self.TaskException(err_msg)
        return True

    def finished(self):
        ''' A method that will check if this Task has finished executing,
            and check to ensure that execution was succesfull. If this Task
            hasn't been started and hasn't been skipped, it will return false.
            Otherwise, this method will query the subproccess created by
            self.start(). If the subprocess is still runnning, this method
            will return False. If the subprocess has set its exit code, the
            exit code will be checked by self.error_check, then this method
            will check to make sure that all targets of this task were
            created. If all checks are passed, self.finished will return True.
            In the event that execution of this task fails for some reason,
            all targets will be renamed to be [t+'.partial' for t in targets]
            to indicate that the output is unlikely to be complete.
        '''
        if(self.soft_finished_status):
            return True
        if(self.process is None):
            return False
        self.process.poll()
        exit_code = self.process.returncode
        if(exit_code is None):
            cur_run_time = float(time.time() - self.start_time) / 60
            if(cur_run_time > self.max_wall_time):
                err_mess = ('Task {0!s} has been running for greater than its maximum wall time, '
                            '{1!s}m, and has been aborted. This is likely an external error and '
                            'trying again is recommended.').format(self.name, self.max_wall_time)
                self.killRun()
                raise self.TaskException(err_mess)
            return False

        else:
            self.exit_code = exit_code
            self.close_files()
            if(self.error_check(exit_code)):
                error_message = ('The task {0!s} seems to have failed with exit_code {1!s}.'
                                 ).format(self.name, exit_code)
                if(self.stderr is not None):
                    error_message += ' Consult {0!s} for more info.'.format(self.stderr)
                inst = self.ExitCodeException(error_message)
                inst.exit_code = exit_code
                self.rename_targets()
                raise inst
            for t in self.targets:
                if(not os.path.exists(t)):
                    error_message = ('Task {0!s} encountered an unexpected error resulting in it'
                                     ' failing to produce its target files.').format(self.name)
                    if(self.stderr is not None):
                        error_message += ' Consult {0!s} for more info.'.format(self.stderr)
                    raise self.TaskException(error_message)
            return True

    def run(self, delay=1):
        ''' A convenience method that alllows for the execution of a task
            serially. This task will be started by a call to self.start(),
            then every delay seconds, a call to self.finished() will check to
            see if the command has finished executing. self.run() will return
            upon completion of the command. 
        '''
        self.start()
        while(not self.finished()):
            time.sleep(delay)

    def killRun(self):
        ''' Calling this method will safely stop the execution of this task.
        '''
        try:
            self.process.kill()
        except:
            pass
        self.close_files()
        self.rename_targets()

    def close_files(self):
        ''' This method is responsible for closing all files opened by a call to self.start.
        '''
        for f in self.opened_files:
            f.close()
        self.opened_files = []

    def rename_targets(self):
        ''' This method will check for each t in self.targets, if t exists.
            If it does, t will be renamed to be t+".partial". This method is
            called whenever this task halts unexpectedly or in error.
        '''
        for f in self.targets:
            if(os.path.exists(f)):
                os.rename(f, f + '.partial')

    def skipable(self):
        ''' Method uses by Supervisor objects to determine if this task needs
            to be executed. This method will return True if this task is
            skippable and so does not need to be executed. A task is skippable
            if and only if all dependencies are satisfied, all dependencies
            have been skippable, and all targets exist.
        '''
        if(self.exit_code is not None):
            return False
        if(self.soft_finished_status):
            return True
        for t in self.dependencies:
            if(isinstance(t, Task) or isinstance(t, Supervisor)):
                if(not t.skipable()):
                    return False
        if(self.targets == []):
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
    ''' Object that will intelligently run Task objects or other Supervisor
        objects.
    '''

    STATE_INITIALIZED = 'initialized'
    STATE_ERR = 'failed'
    STATE_SKIPPED = 'skipped'
    STATE_FINISHED = 'executed'
    STATE_RUNNING = 'started'
    STATE_REMOVED = 'removed'

    def __init__(
      self, tasks=[], dependencies=[], cpu=float('inf'), name='Supervisor',
      delay=1, force_run=False, email=None, email_interval=30, log=None):
        ''' the __init__ for the Superviosr object has nine paramters
            described below.
            tasks -
                A list of Task or Supervisor objects that this Supervisor will
                manage. Default = []
            dependencies -
                A list of dependencies for this Supervisor. This works in the
                same way as the dependencies for a Task object. Default = []
            cpu -
                The cpu cap for this supervisor. When you run this supervisor,
                the supervisor will run the maximimun number of Tasks that it
                can in parrelel such that the sum of Task.cpu for each of those
                tasks is below this supervisor's cpu. Default = float('inf')
            name -
                A handy way to keep track of your supervisor if you are using
                multiple. name will be used during exception help messages and
                when printing the supervisor. Default = 'Supervisor'
            delay -
                While running tasks, this superrvisor will weight delay seconds
                in between each execution cycle. Default = 1
            force_run -
                A flag that determines whether the Supervisor will attempt to
                skip execution of tasks that it thinks are skippable. A task is
                skippable if Task.skippable() returns True. Default = False
            email -
                This supervisor will send emails updating you on execution of
                the pipeline if an email is given. This supervisor will send
                an email when the pipeline finishes, and will send emails at
                most once per email_interval minutes while running.
                Default = None
            email_interval -
                The delay between emails. Sent by the Supervisor. Only relevant
                if email is used.
            log -
                The path to the log file that the Supervisor will generate. The
                log file is plain text and details the execution of each task.
                Default = name+'.run_log'
        '''
        self.cpu = cpu
        self.name = name
        self.delay = delay
        self.dependencies = dependencies
        self.force_run = force_run
        self.email = email
        self.email_interval = email_interval * 60
        self.last_email = time.time()
        self.log_path = log if(log is not None) else name + '.run_log'
        self.log_str = ''
        self.task_map = {}
        self.task_status = {}
        self.errors = []
        self.targets = []
        self.tasks = set()
        for t in tasks:
            self.add_task(t)

    def run(self):
        ''' Use this funciton to get a Supervisor to run all tasks that it is
            managing. Execution occurs inside of a large loop with multiple
            steps. First, this supervisor checks all running tasks to see if
            they are finished. Those tasks that are finished are handled and
            removed from processing Next the Supervisor checks to see if any
            tasks can be started by checking their dependencies and making
            sure executing the task won't result in the SUpervisor exceeding
            the cpu_cap. Then, if force_run is False, the supervisor will
            check to see if the task can be skipped. If the task can be
            skipped, it will be, otherwise the task will be started. The loop
            breaks when their are no tasks left that can be executed and no
            tasks currently being executed.
        '''
        self.log_file = open(self.log_path, 'w', 1)
        self.last_email = time.time()
        self.tasks_to_run = set([t for t in self.tasks])
        self.tasks_running = set()
        cur_cpu = 0
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
                    elif(not self.force_run and t.skipable()):
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
            self.send_email('', subject='MMT Finished')

    def __removeTaskPath__(self, task):
        ''' Helper function that removes tasks that can no longer be executed
            from the sueprvisors list of tasks that need o be executed.
            Basically, anything down stream of task is removed.
        '''
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
                    flag = True

    def log(self, message):
        ''' Convienvince funciton that allows for the SUpervisor to write
            messages to stdout and the supervisors log file
        '''
        self.log_str += message
        self.log_file.write(message)
        if(check_foreground()):
            print(message)
        if(time.time() > self.last_email + self.email_interval):
            self.last_email = time.time()
            self.send_email(self.log_str, 'MMT Running Update')
            self.log_str = ''

    def finished(self):
        ''' A method called from within run. If a task is dependant on a this
            supervisor, then this supervisor will query all of its tasks to
            determine if the supervisor is finished executing. A Supervisor
            is finished if and only if all of the Supervisor's tasks are
            finsihed.
        '''
        for t in self.tasks:
            if(not t.finished()):
                return False
        return True

    def skipable(self):
        ''' A method called from within run. If a task is dependant on a this
            supervisor, then this supervisor will query all of its tasks to
            determine if the supervisor doesn't need to be executed. A
            Supervisor is skippable if and only if all of the Supervisor's
            tasks are skipable.
        '''
        for t in self.tasks:
            if(not t.skipable()):
                return False
        return True

    def send_email(self, message, subject=''):
        ''' A convienince function that allows for the sending of emails to
            the Supervisor's self.email adress. If self.email is none, the
            method does nothing. Otherwise, Supervisor will start a subprocess
            that sends the email with the specified subject and message to
            self.email. subject defaults to ''. The funciton will always return
            None.
        '''
        if(self.email is None):
            return
        else:
            # message = ''.join([c if(c!='\n') else '\\n'for c in message])
            cmd = "echo '{0!s}' | mail -s '{1!s}' '{2!s}'".format(message, subject, self.email)
            subprocess.call(cmd, shell=True)

    def killRun(self):
        ''' safely stops all running tasks and halts the run.
        '''
        for t in self.tasks_running:
            t.killRun()
        self.running = []
        sys.exit(0)

    def add_task(self, task):
        ''' Allows for a task or supervisor to be added to this supervisor.
            A supervisor is added to this supervisor by adding all tasks
            that the supervisor managed.
        '''
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


if(__name__ == '__main__'):
    t1 = Task('python ..\\..\\test.py 4.4', dependencies=[], name='t1')
    t2 = Task('python ..\\..\\test.py 6', dependencies=[t1], name='t2')

    def t2_dep():
        t2.command = 'python ..\\..\\test.py 8'
        return True

    t2.dependencies.append(t2_dep)
    t3 = Task('python ..\\..\\test.py 3', dependencies=[], name='t3')
    t4 = Task('python ..\\..\\test.py 2', dependencies=[t3], name='t4')
    t5 = Task('python ..\\..\\test.py 1', dependencies=[], name='t5')
    t6 = Task('python ..\\..\\test.py 7', dependencies=[], name='t6')
    s1 = Supervisor([t1, t2], name='s1')
    s2 = Supervisor([t3, t4], name='s2', dependencies=[s1])
    s3 = Supervisor([s1, s2], name='s3')
    s4 = Supervisor([t5, t6], name='s4', dependencies=[t1])
    s = Supervisor([s3, s4], force_run=False)
    s.run()

    # help(Task)
    # help(Supervisor)
