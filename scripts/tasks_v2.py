import os
import sys
import time
import signal
import subprocess
import pickle
from functools import reduce
from itertools import chain


class Task:

	class ExitCodeException(Exception):
		exit_code = None

	class TaskException(Exception):
		pass


	def __init__(self,command,dependencies=[],targets=[],cpu=1,name='Anonymous_Task',stderr=None,stdout=None,error_check=None):
		try:
			if(stderr!=None):
				f = open(stderr,'a')
				f.close()
			if(stdout!=None):
				f = open(stdout,'a')
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

	
	def checkDependencies(self):
		for d in self.dependencies:
			if( isinstance(d,Task) ):
				try:
					if( not d.finished()):
						return False
				except Task.ExitCodeException:
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
		temp = subprocess.Popen(self.command,shell=True,stdout=out,stderr=err)
		self.process = temp


	def finished(self):
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


	def __str__(self):
		return self.name

	def __repr__(self):
		return self.__str__()


class Supervisor:

	STATE_INITIALIZED = 'initialized'
	STATE_ERR ='failed'
	STATE_SKIPPED = 'skipped'
	STATE_FINISHED = 'executed'
	STATE_RUNNING = 'started'
	STATE_REMOVED = 'removed'


	def __init__(self,tasks=[],dependencies=[],cpu=12,name='Supervisor',delay=1,force_run=False,email=None,email_interval=30,log=None):
		self.history_log_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'.task_log')
		self.tasks = set(tasks)
		self.running = set([])
		self.cpu = cpu
		self.name = name
		self.delay = delay
		self.dependencies=dependencies
		self.force_run=force_run
		self.email=email 
		self.email_interval=email_interval
		self.last_email = time.time()
		self.log_path = log if(log!=None) else name+'.run_log'
		#self.log_file = open(self.log_path,'w',1)
		self.log_str = ''
		self.task_status = {t:{'state':Supervisor.STATE_INITIALIZED,'exit_code':None,'message':None,'start':None,'stop':None}  for t in tasks}
		self.errors = []
		self.targets = [t for t in chain(*[x.targets for x in tasks])]

	
	def run(self):
		self.log_file = open(self.log_path,'w',1)
		self.last_email=time.time()
		self.checkTasks()
		cur_cpu = 0
		signal.signal(signal.SIGTERM,lambda *x : self.killRun() )
		history_update = set()
		try:
			while(len(self.tasks)>0 or len(self.running)>0):
				for t in [task for task in self.running]:
					try:
						if(t.finished()):
							self.running.remove(t)
							cur_cpu-=t.cpu
							self.task_status[t]['stop'] = int(time.time())
							self.task_status[t]['state'] = Supervisor.STATE_FINISHED
							self.task_status[t]['exit_code'] = t.exit_code
							m,s = divmod(self.task_status[t]['stop']-self.task_status[t]['stop'],60)
							h,m = divmod(m,60)
							self.task_status[t]['message'] = 'Completed in {0!s}h {1!s}m {2!s}s'.format(h,m,s)
							history_update.add(t.command)
							self.log(t.name+':'+self.task_status[t]['state']+':'+time.asctime()+'\n\n')
					except (Task.ExitCodeException,Task.TaskException) as inst:
						self.errors.append(inst)
						self.running.remove(t)
						cur_cpu-=t.cpu
						self.task_status[t]['stop'] = int(time.time())
						self.task_status[t]['state'] = Supervisor.STATE_ERR
						self.task_status[t]['exit_code'] = t.exit_code
						m,s = divmod(self.task_status[t]['stop']-self.task_status[t]['start'],60)
						h,m = divmod(m,60)
						self.task_status[t]['message'] = 'Failed in {0!s}h {1!s}m {2!s}s'.format(h,m,s)
						self.log(t.name+':'+self.task_status[t]['state']+':'+time.asctime()+'\n\n')						
						self.__removeTaskPath__(t)
				
				for t in [task for task in self.tasks]:
					if(t.cpu+cur_cpu<=self.cpu and t.checkDependencies()):
						cur_cpu+=t.cpu
						self.running.add(t)
						self.tasks.remove(t)
						t.start()
						self.task_status[t]['state'] = Supervisor.STATE_RUNNING
						self.task_status[t]['start'] = int(time.time())
						self.log(t.name+':'+self.task_status[t]['state']+':'+time.asctime()+'\n\n')
						self.log_file.write(t.command+'\n\n')

				if(len(self.running)==0 and len(self.tasks)!=0):
					break
				time.sleep(self.delay)
			if(self.errors!=[] or (len(self.running)==0 and len(self.tasks)!=0)):
				err_str = '\n\n'
				if(len(self.running)==0 and len(self.tasks)!=0):
					err_str+= 'Unable to resolve dependencies during execution of '+self.name
					err_str+='. The following tasks could not be executed:\n'
					err_str+= '\n'.join(['\t'+t.name+':dependencies-'+str([d.name for d in t.dependencies]) for t in self.tasks])
				if(self.errors!=[]):
					err_str += '\nEncountered an unexpected Error in the following tasks:\n'
					for t in self.task_status:
						if(self.task_status[t]['state']==Supervisor.STATE_ERR):
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
			self.send_email('',subject='Pipeline Finished')


	def __removeTaskPath__(self,task):
		flag = True
		removed = set([task])
		check_intersection = lambda s,l : any(e in s for e in l)
		while(flag):
			flag = False
			for t in [t for t in self.tasks]:
				if(check_intersection(removed,t.dependencies)):
					self.tasks.remove(t)
					removed.add(t)
					self.task_status[t]['message'] = 'Never Started'
					self.task_status[t]['state'] = Supervisor.STATE_REMOVED
					flag=True

	
	def checkTasks(self):
		self.getTasks()		
		self.task_status = {t:{'state':Supervisor.STATE_INITIALIZED,'exit_code':None,'message':None,'start':None,'stop':None}  for t in self.tasks}
		for t in self.tasks:
			if(t.cpu>self.cpu):
				error_message = 'Task '+t.name+' has a higher cpu than this supervisor, '+self.name+'.\n'
				error_message+= 'You must increase this supervisor\'s cpu or decrease '+t.name+'\'s cpu.'
				raise Exception(error_message)
		if(self.force_run):
			return
		r_graph = {t:set([d for d in t.dependencies if(isinstance(d,Task) or isinstance(d,Supervisor))]) for t in self.tasks}
		f_graph = {t:set([d for d in r_graph if(t in r_graph[d]) ]) for t in self.tasks}
		task_required = {t:False for t in self.tasks}
		task_history = self.get_history()
		flag = True
		check_targets = lambda x : not all(os.path.exists(z) for z in x)
		while(flag):
			flag = False
			for t in self.tasks:
				if(task_required[t]):
					continue
				if(t.command not in task_history or check_targets(t.targets)):
					task_required[t] = True
					flag = True
				for d in r_graph[t]:
					if(task_required[d]):
						task_required[t] = True
						flag = True
		for t in task_required:
			if(not task_required[t]):
				self.tasks.remove(t)
				self.task_status[t]['state'] = Supervisor.STATE_SKIPPED
				self.log(t.name+':'+self.task_status[t]['state']+'\n')
				for d in f_graph[t]:
					if(t in d.dependencies):
						d.dependencies.remove(t)


	def getTasks(self):
		r_graph = {t:set([d for d in t.dependencies if(isinstance(d,Task) or isinstance(d,Supervisor))]) for t in self.tasks}
		f_graph = {t:set([d for d in r_graph if(t in r_graph[d]) ]) for t in self.tasks}
		flag = True
		while(flag):
			flag=False
			for task in [t for t in self.tasks]:
				if(isinstance(task,Supervisor)):
					flag=True
					sub_tasks = task.tasks
					for elem in sub_tasks:
						elem.dependencies.extend(task.dependencies)
						self.tasks.add(elem)
					for dependant_task in f_graph[task]:
						dependant_task.dependencies.remove(task)
						dependant_task.dependencies.extend(sub_tasks)
					self.tasks.remove(task)
				r_graph = {h:set([d for d in h.dependencies if(isinstance(d,Task) or isinstance(d,Supervisor))]) for h in self.tasks}
				f_graph = {h:set([d for d in r_graph if(h in r_graph[d]) ]) for h in self.tasks}
		return self.tasks


	def get_history(self):
		try:
			history_file = open(self.history_log_path,'rb')
			ret = pickle.load(history_file)
			history_file.close()
		except (KeyboardInterrupt, SystemExit):
			raise
		except:
			ret = set()
		return ret


	def write_history(self,update):
		history = self.get_history()
		history = history.union(update)
		history_file = open(self.history_log_path,'wb')
		pickle.dump(history,history_file)
		history_file.close()


	def log(self,message):
		self.log_str+=message
		self.log_file.write(message)
		print(message)
		if(time.time() > self.last_email+self.email_interval):
			self.last_email = time.time()
			#self.send_email(self.log_str,'Pipeline Running Update')
			self.log_str=''


	def send_email(self,message,subject=''):
		if(self.email==None):
			return
		else:
			message = ''.join([c if(c!='\n') else '\\n'for c in message])
			cmd = "echo '{0!s}' | mail -s '{1!s}' '{2!s}'".format(message,subject,self.email)
			subprocess.call(cmd)


	def killRun(self):
		for t in self.running:
			t.killRun()
		self.running = []


	def __str__(self):
		return self.name

	def __repr__(self):
		return self.__str__()



if(__name__=='__main__'):
	t1 = Task('py3 ..\\..\\test.py 4',dependencies=[],name='t1')
	t2 = Task('py3 ..\\..\\test.py 6',dependencies=[t1],name='t2')
	t3 = Task('py3 ..\\..\\test.py 3',dependencies=[],name='t3')
	t4 = Task('py3 ..\\..\\test.py 2',dependencies=[t3],name='t4')
	t5 = Task('py3 ..\\..\\test.py 1',dependencies=[],name='t5')
	t6 = Task('py3 ..\\..\\test.py 7',dependencies=[],name='t6')
	s1 = Supervisor([t1,t2],name='s1')
	s2 = Supervisor([t3,t4],name='s2',dependencies=[s1])
	s3 = Supervisor([s1,s2],name='s3')
	s4 = Supervisor([t5,t6],name='s4',dependencies=[t1])
	s = Supervisor([s1,s3,s4],force_run=True)
	s.run()


