#!/usr/bin/python

#---AUTOMACS (AMX) "TOOLS" MODULE
#---created 2014.08.13 by ryan bradley

"""
The tools module provides supporting functions for the other modules in order to simplify the code. These 
consist primarily of functions which interact with the filesystem either via Python libraries such as 
``shutil`` or via Bash using the ``subprocess`` or ``sys`` modules. Modules import necessary functions with 
the following code. ::

	from amx.tools import call,tee,checkout
"""

import time
import datetime
import subprocess
import sys,os
import glob,re
import errno
import shutil

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

class tee(object):
	'''
	Routes print statements to the screen and a log file.
	
	Routes output to multiple streams, namely a log file and stdout, emulating the linux "tee" function. 
	Whenever a new log file is desired, use the following code to replace ``sys.stdout`` and route all print
	statements to both the screen and the log file. ::
		
		sys.stdout = tee(open(self.rootdir+'log-script-master','a',1))
		
	Initialize the object with a file handle for the new log file. It is possible to run ``tee`` multiple
	times in order to redirect ``print`` output to a new file. The new object checks to see if 
	``sys.stdout`` is a tee object or the "real" stream, and rolls both into the new object.
	'''
	def __init__(self, *files,**kwargs):
		#---if sys.stdout is already a tee object, then just steal its stdout member
		if str(sys.stdout.__class__) == "<class 'amx.tools.tee'>": self.stdout = sys.stdout.stdout
		#---otherwise set stdout from scratch
		else: 
			if 'error' in kwargs.keys() and kwargs['error'] == True: self.stdout = sys.stderr
			else: self.stdout = sys.stderr
		self.files = files
		if 'error' in kwargs.keys() and kwargs['error'] == True: self.error = True
		else: self.error = False
	def write(self, obj): 
		'''
		The write function here emulates the write functions for both files and the standard output stream
		so that the tee object will always write to both places.
		'''
		st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d.%H%M') if not self.error else ''
		if obj != '\n': self.stdout.write(st+' ')
		self.stdout.write(obj)
		for f in self.files:
			if obj != '\n': f.write(st+' ')
			f.write(obj)

def call(command,logfile=None,cwd=None,silent=False,inpipe=None):
	'''
	Wrapper for system calls in a different directory with a dedicated log file.
	'''
	if inpipe != None:
		output = open(('' if cwd == None else cwd)+logfile,'wb')
		if type(command) == list: command = ' '.join(command)
		p = subprocess.Popen(command,stdout=output,stdin=subprocess.PIPE,stderr=output,cwd=cwd,shell=True)
		catch = p.communicate(input=inpipe)[0]
	else:
		if type(command) == list: command = ' '.join(command)
		if logfile != None:
			output = open(('' if cwd == None else cwd)+logfile,'wb')
			if type(command) == list: command = ' '.join(command)
			if not silent: print 'executing command: "'+str(command)+'" logfile = '+logfile
			try:
				subprocess.check_call(command,
					shell=True,
					stdout=output,
					stderr=output,
					cwd=cwd)
			except: 
				if logfile[-3:] == '-cg' and re.search('mdrun-em',logfile):
					print 'warning: failed conjugate gradient descent but will procede'
				else: raise Exception('except: execution error')
			output.close()
		else: 
			if not silent: print 'executing command: "'+str(command)+'"'
			if str(sys.stdout.__class__) == "<class 'amx.tools.tee'>": stderr = sys.stdout.files[0]
			else: stderr = sys.stdout
			try: subprocess.check_call(command,shell=True,stderr=stderr,cwd=cwd)
			except:
				if logfile[-3:] == '-cg' and re.search('mdrun-em',logfile):
					print 'warning: failed conjugate gradient descent but will procede'
				else: raise Exception('except: execution error')

def checkout(command,cwd=None):
	'''
	Return the result of a shell command.
	
	This command emulates the ``subprocess.check_output`` function for systems with older versions of that 
	module. If ``subprocess.check_output`` is unavailable, it routes the ``STDOUT`` from the command to
	a temporary text file and then reads that.
	'''
	try:
		result = subprocess.check_output(' '.join(command),shell=True,cwd=cwd)
	except AttributeError:
		if type(command) != list: raise Exception('except: command must be a list for old-school Popen')
		if command[0] == 'awk': command[1] = re.sub("'","",command[1])
		result = subprocess.Popen(command,stdout=subprocess.PIPE,cwd=cwd).communicate()[0].strip('\n')
	return result

def copy(source,destination):
	'''
	Generic copy function for handling files or directories.
	
	This copy function can copy directories as long as the destination folder is fully specified. It will
	also copy a folder's contents to another folder directly if the Bash-style asterisk wildcard is found
	in the source argument. It also copies individual files.
	'''
	#---check if wildcard in the source in which case use glob
	if '*' in source:
		if 0: print glob.glob(source)
		for filename in glob.glob(source):
			shutil.copy(filename,destination)
	#---otherwise try a copytree and if that fails just use copy
	else:
		try: shutil.copytree(source,destination)
		except OSError as exc:
			if exc.errno == errno.ENOTDIR: shutil.copy(source,destination)
			else: raise

def LastFrame(prefix,rootdir,gmxpaths):
	'''
	This function is used to retreive the last frame of a particular simulation.\n
	The user supplies a filename prefix and root directory (via kwargs `prefix` and `rootdir`) and this
	function will check for the corresponding <prefix>.gro file. If it's absent, it checks for a <prefix>.xtc 
	file and extract the last frame from it.
	'''
	if not os.path.isfile(rootdir+'/'+prefix+'.gro'):
		print 'configuration is missing so we will extract the last frame'
		cmd = [gmxpaths['gmxcheck'],
			'-f '+prefix+'.xtc']
		call(cmd,logfile='log-gmxcheck-'+prefix,cwd=rootdir)
		lastframe = int(checkout(["awk","'/Step / {print $2}'",
			'log-gmxcheck-'+prefix],cwd=rootdir).strip())
		ts = float(checkout(["awk","'/Step / {print $3}'",
			'log-gmxcheck-'+prefix],cwd=rootdir).strip())
		#---note that rpb added a rounding function after an error in an AAMD test 2014.08.24
		lasttime = float(int(float(lastframe)*ts))
		lasttime = round(lasttime/10)*10
		print 'last viable frame was at '+str(lasttime)+' ps'
		print 'retrieving that frame'
		cmd = [gmxpaths['trjconv'],
			'-f '+prefix+'.xtc',
			'-o '+prefix+'.gro',
			'-s '+prefix+'.tpr',
			'-b '+str(lasttime),
			'-e '+str(lasttime),]
		call(cmd,logfile='log-trjconv-'+prefix,cwd=rootdir,inpipe='0\n')
	else: print 'configuration file is already available'

