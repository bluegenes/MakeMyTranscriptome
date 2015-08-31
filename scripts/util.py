'''	Author : Nolan Hartwick
	Date : 6/5/15
	Overview : Module contains random useful funcitons.
	Contents :
		dp : Used to convert recursive function to dynamic function.
		logger : Adds logging funcitonality to a function.
		memoize : Used to convert a function to a memoized function.
		Wrapper : A decorator for decorators that fixes doc strings / names.
'''


from sys import stdout
from time import time


def wrapper(g,dec_name):
	'''	A parameterized decorator function designed to be used on the 
		functions inside and returned by decorators. Wrapper should always 
		be used with parenthesis. It will set a decorated functions name to 
		be g.name + _ + dec_name. It will also set a decorated functions doc 
		string to be equal to g.doc + newline + 'Decorated with ' + dec_name.
	'''
	def decorator(f):
		def ret(*args,**kargs):
			return f(*args,**kargs)
		try:
			ret.__name__ = g.__name__ + '_' + dec_name
			ret.__doc__ = g.__doc__ + '\nDecorated with ' + dec_name
		except (KeyboardInterrupt, SystemExit):
			raise
		except:
			pass
		ret.original_function = g
		return ret
	return decorator


def dp(f):
	'''	Decorator that will dynamically program a function f. To access the 
		table used to store the call and return values, simply set a variable 
		to be equal to f.table prior to executing f. After executing f, your 
		variable will now provide access to the filled out table. A subsequent 
		call to f will not modify the variable or the table associated with it.
		The table is a dict of key value pairs in which the keys are a set of 
		args passed to f, either on initial call or in any recursive call, and 
		the values are f(args).
	'''
	@wrapper(f,'dp')
	def ret(*args,**kargs):
		flag = ret._top
		if(flag):
			ret._top = False
		if( args not in ret.table):
			ret.table[args] = f(*args, **kargs)
		r = ret.table[args]
		if(flag):
			ret.table = {}
			ret._top = True
		return r
	ret.table = {}
	ret._top = True
	return ret


def memoize(f):
	'''	Decorator that will convert a function f into a memoized function. To 
		access the call table, simply use f.table. 
	'''
	@wrapper(f,'memoize')
	def ret(*args,**kargs):
		if(args not in ret.table):
			ret.table[args] = f(*args,**kargs)
		return ret.table[args]
	ret.table = {}
	return ret


def logger(seperator=' : ',outfile=stdout):
	'''	A parameterized decorator function for logging. When used as a 
		decorator, always use parenthesis. The decorated function will now 
		write its name, positional arguments, keyword arguments, and return 
		values seperated by seperator to the opened writeable file outfile 
		when the function f returns. If f raises an error, everything will 
		be printed as normal with the exception of the return value which 
		will be replaced with 'Raised Error' 
	'''
	def dec(f):
		@wrapper(f,'logger')
		def ret(*args, **kargs):
			p_string = f.__name__+seperator+str(args)+seperator+str(kargs)+seperator
			try:
				ret = f(*args,**kargs)
				outfile.write(p_string+str(ret)+'\n')
				return ret
			except:
				outfile.write(p_string+'Raised Error\n')
				raise
		return ret
	return dec


if(__name__=='__main__'):
	
	
	#@logger()
	def fib(n):
		'''	 This is Fib.
		'''
		if(n<=2):
			return 1
		return fib(n-1)+fib(n-2)


	def gloabal_allign(i,j,seq1='',seq2='',score=lambda x,y : 1 if(x==y) else -1, indel=-1):
		if(i==0 and j==0):
			return 0
		elems = []
		if(j>0):
			elems.append(gloabal_allign( i, j-1, seq1=seq1, seq2=seq2, score=score, indel=indel) + indel)
		if(i>0):
			elems.append(gloabal_allign( i-1, j, seq1=seq1, seq2=seq2, score=score, indel=indel) + indel)
		if(i>0 and j>0):
			elems.append(gloabal_allign( i-1, j-1, seq1=seq1, seq2=seq2, score=score, indel=indel) + score( seq1[i-1], seq2[j-1]))
		return max(elems)
	
	start = time()
	fib(37)
	print(time()-start)
