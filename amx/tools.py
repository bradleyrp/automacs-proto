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
import socket
from copy import deepcopy
if os.path.isfile('settings'): execfile('settings')
elif os.path.isfile('../settings'): execfile('../settings')
if 'USER' in os.environ: username = os.environ['USER']
else: username = re.findall('/home/(\w+)',os.environ['HOME']).pop()
extrasets = os.path.expanduser('/home/%s/.automacs.py'%username)
if os.path.isfile(extrasets): execfile(extrasets)

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

#---classic python argsort
def argsort(seq): return [x for x,y in sorted(enumerate(seq), key = lambda x: x[1])]

class tee(object):

	"""
	Routes print statements to the screen and a log file.
	
	Routes output to multiple streams, namely a log file and stdout, emulating the linux "tee" function. 
	Whenever a new log file is desired, use the following code to replace ``sys.stdout`` and route all print
	statements to both the screen and the log file. ::
		
		sys.stdout = tee(open(self.rootdir+'log-script-master','a',1))
		
	Initialize the object with a file handle for the new log file. It is possible to run ``tee`` multiple
	times in order to redirect ``print`` output to a new file. The new object checks to see if 
	``sys.stdout`` is a tee object or the "real" stream, and rolls both into the new object.
	"""

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
		"""
		The write function here emulates the write functions for both files and the standard output stream
		so that the tee object will always write to both places.
		"""
		st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d.%H%M') if not self.error else ''
		if obj != '\n': self.stdout.write(st+' ')
		self.stdout.write(obj)
		for f in self.files:
			if obj != '\n': f.write(st+' ')
			f.write(obj)

def call(command,logfile=None,cwd=None,silent=False,inpipe=None,suppress_stdout=False):

	"""
	Wrapper for system calls in a different directory with a dedicated log file.
	"""

	if inpipe != None:
		if logfile == None: output=None
		else: output = open(('' if cwd == None else cwd)+logfile,'wb')
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
				else: raise Exception('except: BASH execution error (probably GROMACS)!\nsee '+cwd+logfile)
			output.close()
		else: 
			if not silent: print 'executing command: "'+str(command)+'"'
			if str(sys.stdout.__class__) == "<class 'amx.tools.tee'>": stderr = sys.stdout.files[0]
			else: stderr = sys.stdout
			try: 
				if suppress_stdout: 
					devnull = open('/dev/null','w')
					subprocess.check_call(command,shell=True,stderr=devnull,cwd=cwd,stdout=devnull)
				else: subprocess.check_call(command,shell=True,stderr=None,cwd=cwd)
			except: 
				raise Exception('except: BASH execution error (probably GROMACS)!\ncommand: '+\
					command+'\ncwd: '+cwd)

def checkout(command,cwd=None):

	"""
	Return the result of a shell command.
	
	This command emulates the ``subprocess.check_output`` function for systems with older versions of that 
	module. If ``subprocess.check_output`` is unavailable, it routes the ``STDOUT`` from the command to
	a temporary text file and then reads that.
	"""

	try:
		result = subprocess.check_output(' '.join(command),shell=True,cwd=cwd)
	except AttributeError:
		if type(command) != list: raise Exception('except: command must be a list for old-school Popen')
		if command[0] == 'awk': command[1] = re.sub("'","",command[1])
		result = subprocess.Popen(command,stdout=subprocess.PIPE,cwd=cwd).communicate()[0].strip('\n')
	return result

def copy(source,destination):

	"""
	Generic copy function for handling files or directories.
	
	This copy function can copy directories as long as the destination folder is fully specified. It will
	also copy a folder's contents to another folder directly if the Bash-style asterisk wildcard is found
	in the source argument. It also copies individual files.
	"""

	#---check if wildcard in the source in which case use glob
	if '*' in source:
		for filename in glob.glob(source): 
			if os.path.isdir(filename): shutil.copytree(filename,destination+os.path.basename(filename))
			else: shutil.copy(filename,destination)
	#---otherwise try a copytree and if that fails just use copy
	else:
		try: shutil.copytree(source,destination)
		except OSError as exc:
			if exc.errno == errno.ENOTDIR: shutil.copy(source,destination)
			else: raise

def multiresub(switches,text):

	"""
	Function which performs a multiple regex substitution on some text.\n
	Primarily used in controller to prepare ad hoc bash scripts.
	"""

	regex = re.compile("(%s)" % "|".join(map(re.escape,switches.keys())))
	return regex.sub(lambda mo: switches[mo.string[mo.start():mo.end()]],text) 
	
def confirm():

	"""
	Generic function to check with the user.
	"""
	
	go = True if raw_input("%s (y/N) " % 'continue?').lower() == 'y' else False
	if not go:
		print 'aborting' 
		return False
	sure = True if raw_input("%s (y/N) " % 'confirmed?').lower() == 'y' else False
	if go and sure: return True
	
def status(string,start=None,i=None,looplen=None,blocked=False):
	
	"""
	Print status to the screen also allows for re-writing the line. Duplicated in the membrain library.
	"""
	
	#---note: still looking for a way to use the carriage return for dynamic counter without
	#---...having many newlines printed to the file. it seems impossible to use the '\r' + flush()
	#---...method with both screen and file output, since I can't stop the buffer from being written
	#---the blocked flag will remove tabs so that large multiline text doesn't awkwardly wrap in the code
	if blocked: string = niceds(text)
	#---display a refreshable string
	if start == None and looplen == None and i != None:		
		print '\r'+string+'  ...  '+str(i+1).rjust(7)+'/'+str(looplen).ljust(8)+'\n',
	elif start == None and looplen != None and i != None:		
		if i+1 == looplen:
			print '\r'+string+'  ...  '+str(i+1).rjust(7)+'/'+str(looplen).ljust(8)+'\n',
		#---if a logfile has been defined, this output is destined for a file in which case suppress counts
		elif i+1 != looplen and ('logfile' not in globals() or logfile == None):
			print '\r'+string+'  ...  '+str(i+1).rjust(7)+'/'+str(looplen).ljust(8),
			sys.stdout.flush()
	#---estimate the remaining time given a start time, loop length, and iterator
	elif start != None and i != None and looplen != None and ('logfile' not in globals() or logfile == None):
		esttime = (time.time()-start)/(float(i+1)/looplen)
		print '\r'+string.ljust(20)+str(abs(round((esttime-(time.time()-start))/60.,1))).ljust(10)+\
			'minutes remain',
		sys.stdout.flush()
	#---standard output
	else: print string
	
def chain_steps():
	
	"""
	Check for previous steps to enable chaining.
	"""
	
	#---intervene here to check if this is an add-on singleton script
	#---...for example, if you want a restart of a previous procedure
	#---...this section allos you to chain procedures together
	#---...note that no backup of previous scripts is made here so the script-master is overwritten
	for root,dirnames,filenames in os.walk('./'): break
	stepdirs = [i for i in dirnames if re.match('^[a-z][0-9]{1,2}-',i)]
	stepnums = [int(re.findall('^[a-z]([0-9]{1,2})-(.+)',i)[0][0]) for i in stepdirs]
	oldsteps = [stepdirs[i] for i in argsort(
		[int(re.findall('^[a-z]([0-9]{1,2})-(.+)',i)[0][0]) for i in stepdirs])]
	if oldsteps != []:
		if 0: print '\tchaining new steps after old ones'
		if 0: print '\tprevious sequence ended with '+oldsteps[-1]
		startstep = int(re.findall('^[a-z]([0-9]{1,2})-(.+)',oldsteps[-1])[0][0])
	else: startstep = 0
	return startstep,oldsteps
	
def write_steps_to_bash(steps,startstep,oldsteps,extras=None):
	
	"""
	Given a file pointer and a list of steps from the script_dict, this function will write the relevant 
	variables to bash variables. Note that we use this function in script_maker and prep_scripts (which
	is why we place it in a central location here). The steps dictionary does two things: it defines the 
	folders for steps which the execute and preparation scripts require, and also handles special functions
	e.g. looking for previous step names.
	"""

	#---???
	#---only rescripts are in here !!!

	"""
	regex='([a-z,A-Z,0-9_]+)=(.+)';ans = [re.findall(regex,l)[0] for l in lines if re.match(regex,l)];print '\n'.join([str(i) for i in ans])
	"""

	overrides = {}
	regex = '([a-z,A-Z,0-9_]+)=([^$]+)'
	if 'rescript_variables_from' in extras:
		with open('script-'+extras['rescript_variables_from'],'r') as fp: lines = fp.readlines()
		overrides = dict([[i.strip('\n') for i in re.findall(regex,l)[0]]
			for l in lines if re.match(regex,l)])

	bashheader = ''
	#---write variables in the steps entry in the simulation dictionary so bash can see them
	for step in steps.keys(): 
		if step in overrides.keys():
			bashheader += step+'='+overrides[step]+'\n'
			print "OVERRIDING"
			print overrides[step]
		#---intervene to see if this is a step or not
		elif steps[step] != None and re.match('^[a-z]([0-9]{1,2})-(.+)',steps[step]) \
			and ('rescript' not in extras.keys() or not extras['rescript']):
			si,stepname = re.findall('^[a-z]([0-9]{1,2})-(.+)',steps[step])[0]
			#---maintain first letter and advance the step number to append the new step
			newname = steps[step][0]+('%02d'%(int(si)+startstep))+'-'+stepname
			bashheader += step+'='+newname+'\n'
		#---intervene to point to the previous step if the detect_previous_step key is present
		#---...and oldsteps is not empty
		elif step == 'detect_previous_step' and oldsteps != []:
			bashheader += 'detect_previous_step='+oldsteps[-1]+'\n'
		elif step not in ['detect_previous_step']: bashheader += step+'='+steps[step]+'\n'
	if extras != None:
		for key in extras.keys():
			if key not in ['carefultime','rescript']:
				bashheader += key+'='+str(extras[key])+'\n'
	return bashheader
	
def script_maker(target,script_dict,module_commands=None,sim_only=False,stepcount=None,extras=None):

	"""
	Prepare the master script for a particular procedure.
	"""
	if target not in script_dict.keys(): raise Exception('except: unclear make target')
	steps = script_dict[target]['steps']
	script = '\n#---definitions\n'
	startstep,oldsteps = chain_steps()
	#---for time-sensitive cluster execution allow stepcount to truncate the number of steps in sequence
	if 'carefultime' in extras.keys() and extras['carefultime']: stepcount = -1
	#---if boolean the sim_only flag uses only the last, presumably "continuation" segment of the sequence
	#---if list then the sim_only flag supplies the restart script
	script += write_steps_to_bash(steps,startstep,oldsteps,extras=extras)+'\n'
	if type(sim_only) == str: seq_segments = [str(sim_only)]
	elif type(sim_only) == bool and sim_only: seq_segments = script_dict[target]['sequence'][-1:]
	elif stepcount != None: seq_segments = script_dict[target]['sequence'][:stepcount]
	else: seq_segments = script_dict[target]['sequence']
	for segment in seq_segments: 
		#---we add module_commands before each segment to ensure GROMACS paths are set
		if module_commands != None: script += '\n'+module_commands+'\n\n'
		#---if sim_only, then the script is run from within the step folder in which case we remove cd
		if type(sim_only) == str or (type(sim_only)==bool and sim_only):
			segtrim,seglines = [],segment.split('\n')
			for line in seglines:
				if not re.match('^cd\s',line):
					script += line+'\n'
		else: script += segment
	return script
	
def prep_scripts(target,script_dict,extras=None,extra_settings=None,single=None):

	"""
	Execute temporary bash scripts to setup the directories.
	"""
	
	#---extra settings can come from script_dict steps or override via flags passed to make
	#---note that the canonical sets_flags must match those in controller.script
	exset = {} if extra_settings == None else extra_settings
	sets_flags = ['gpu_flag']
	for f in sets_flags:
		if f in script_dict[target]['steps'] and f not in exset:
			exset[f] = script_dict[target]['steps'][f]
	
	startstep,oldsteps = chain_steps()
	if target not in script_dict.keys(): raise Exception('except: unclear make target')
	steps = script_dict[target]['steps']
	for segment in script_dict[target]['prep']:
		fp = open('script-temporary-prep.sh','w')
		fp.write('#!/bin/bash\n\n')
		bashheader = write_steps_to_bash(steps,startstep,oldsteps,extras=extras)		
		fp.write(bashheader)
		fp.write(segment)
		fp.close()
		os.system('bash script-temporary-prep.sh')
		#---any extra_settings are appended to settings.sh files if they appear in a folder 
		#---...labelled by step_name in the steps key of the script_dict
		if extra_settings != None:
			#---extract any step folders that might have settings.sh files
			step_folders = [re.findall('^.+=(.+)$',s)[0] 
				for s in bashheader.split('\n') if any([re.match(r,s) 
				for r in ['^step_','detect_previous_step']])]
			for s in step_folders:
				if os.path.isfile(s+'/settings.sh'):
					with open(s+'/settings.sh','a') as fp: 
						for ex in exset: fp.write(ex+'='+exset[ex]+'\n')
		os.system('rm script-temporary-prep.sh')

def niceblock(text):
	
	"""
	Remove tabs so that large multiline text doesn't awkwardly wrap in the code.
	"""
	
	r"""
	note the previous version of this function returned only the following:
	re.sub('\n([\t])+',(' ' if not newlines else '\n'),re.sub('^\n([\t])+','',text))
	where newlines==True allowed you to keep any explicit newlines
	this function basically just removed all indents so that you could embed blocks 
	of text at a nice tab level and have them print without tabs
	the updated code below removes as many tabs as the first non-newline 
	which makes it suitable for storing code nuggets for printing and execution
	"""
	
	leading_tabs = len(re.findall('\t+(?![\t])',text.strip('\n'))[0])
	tabchop = '\t{%d}(?![\t])'%leading_tabs
	return '\n'.join([(
		re.findall('\t{%d}(.+)'%leading_tabs,line)[0]
		if leading_tabs > 0 and re.match('\t{%d}(.+)'%leading_tabs,line)
		else line
		) for line in text.split('\n')])
		
	
def wordwrap(text,width=60,uniform=False):

	"""
	Basic word wrapper that preserves newlines. Try it with niceblock.
	"""

	#---uniform allows the user to ignore newlines in the text
	if uniform: text = re.sub('\n',' ',text.strip('\n'))
	return reduce(lambda line,word,width=40: '%s%s%s'%(line,
		' \n'[(len(line)-line.rfind('\n')-1+len(word.split('\n',1)[0])>=width)],
		word),text.split(' '))

def noteblock(msg,category='NOTE'):

	"""
	Standard method for printing notes.
	"""

	return '\n'+'\n'.join(['['+category+'] %s'%s for s in 
		wordwrap(niceblock(msg),width=50,uniform=True).split('\n')])+'\n'
	
def latestcheck(last):
	
	"""
	Return the most recept tpr/cpt files in a particular folder.
	"""
	
	for root,dirnames,filenames in os.walk(last): break
	tprs = [i for i in filenames if re.match('^md\.part[0-9]{4}\.tpr$',i)]
	cpts = [i for i in filenames if re.match('^md\.part[0-9]{4}\.cpt$',i)]
	if tprs == [] or cpts == []: return []
	tnum = argsort([int(re.findall('^md\.part([0-9]{4})\.tpr$',i)[0]) for i in tprs])[-1]
	cnum = argsort([int(re.findall('^md\.part([0-9]{4})\.cpt$',i)[0]) for i in cpts])[-1]
	if tprs[tnum][:-4] != cpts[cnum][:-4]: raise Exception('latest cpt/tpr files have different indices')
	return [tprs[tnum],cpts[cnum]]
	
def lastframe(prefix,rootdir,gmxpaths):

	"""
	This function is used to retreive the last frame of a particular simulation.\n
	The user supplies a filename prefix and root directory (via kwargs `prefix` and `rootdir`) and this
	function will check for the corresponding <prefix>.gro file. If it's absent, it checks for a <prefix>.xtc 
	file and extract the last frame from it.
	"""

	if not os.path.isfile(rootdir+'/'+prefix+'.gro'):
		print 'configuration is missing so we will extract the last frame'
		cmd = [gmxpaths['gmxcheck'],
			'-f '+prefix+'.xtc']
		call(cmd,logfile='log-gmxcheck-'+prefix,cwd=rootdir+'/')
		lastframe = int(checkout(["awk","'/Step / {print $2}'",
			'log-gmxcheck-'+prefix],cwd=rootdir+'/').strip())
		ts = float(checkout(["awk","'/Step / {print $3}'",
			'log-gmxcheck-'+prefix],cwd=rootdir+'/').strip())
		t0 = float(checkout(["awk","'/Reading\s+frame\s+0\s+time/ {print $5}'",
			'log-gmxcheck-'+prefix],cwd=rootdir+'/').strip())
		#---note that rpb added a rounding function after an error in an AAMD test 2014.08.24
		lasttime = float(int((float(lastframe)-1)*ts))
		lasttime = round(lasttime/10)*10
		print 'last viable frame was at '+str(lasttime)+' ps'
		print 'retrieving that frame'
		cmd = [gmxpaths['trjconv'],
			'-f '+prefix+'.xtc',
			'-o '+prefix+'.gro',
			'-s '+prefix+'.tpr',
			'-b '+str(lasttime),
			'-e '+str(lasttime),
			'-pbc mol',]
		print ' '.join(cmd)
		call(cmd,logfile='log-trjconv-'+prefix,cwd=rootdir+'/',inpipe='0\n')
	else: print 'configuration file is already available'
	
def get_proc_settings():

	"""
	Determine processor settings from hostname.
	"""

	#---determine location
	hostname = None
	for key in valid_hostnames.keys():
		subproc_failed = False
		try: check_host = re.match(key,subprocess.check_output(['echo $HOSTNAME'],
			shell=True).strip('\n'))
		except: check_host = False
		if check_host or re.match(key,socket.gethostname(),flags=re.I):
			hostname = key
			break
	print '\thostname = '+str(hostname)
	#---if no hostname matches use the local steps
	if hostname == None: system_id = 'local'
	#---check for local gromacs version and update the overrides if necessary
	gmxversion = 5
	try: call('gmx',suppress_stdout=True,silent=True)
	except: gmxversion = 4
	#---check for multiple architectures but select the 
	if hostname in valid_hostnames.keys() and valid_hostnames[hostname] != None and \
		type(valid_hostnames[hostname])==list: 
		arch = valid_hostnames[hostname][0]
	elif hostname in valid_hostnames.keys() and valid_hostnames[hostname] != None:
		arch = valid_hostnames[hostname]
	else: arch = None
	if hostname != None: system_id = hostname+('' if arch == None else '_'+arch)
	if hostname != None and system_id in default_proc_specs.keys():
		print '\thost/architecture settings: '
		for key in default_proc_specs[system_id]: 
			print '\t\t'+key+' = '+str(default_proc_specs[system_id][key])
		proc_settings = default_proc_specs[system_id]
	else: proc_settings = None
	#---in case gromacs 5 is not loaded before the check above we check for modules that match gromacs-5
	if proc_settings != None and \
		'module' in proc_settings and \
		proc_settings['module'] != None and \
		re.match('.+gromacs-5',proc_settings['module']): 
		gmxversion = 5
	#---write gmxpaths.conf
	#---note that you can only change the NPROCS commands passed downstream here
	#---...this means that ./controller make must be rerun to switch processor counts
	fp = open('gmxpaths.conf','w')
	for key in standard_gromacs_commands:
		if gmxversion==5 or (system_id in gmx_overrides.keys() and key in gmx_overrides[system_id].keys()):
			if system_id in gmx_overrides and key in gmx_overrides[system_id]: 
				command_syntax = gmx_overrides[system_id][key]
			else: command_syntax = key
			if proc_settings != None and 'nodes' in proc_settings:
				command_syntax = re.sub('NPROCS',
					str(proc_settings['nodes']*proc_settings['ppn']),command_syntax)
			if gmxversion==5: command_syntax = re.sub(key,gmx_syntax5[key],command_syntax)
			fp.write(key+' '+command_syntax+'\n')
		elif system_id in gmx_suffixes.keys(): 
			fp.write(key+' '+key+gmx_suffixes[system_id]+'\n')
		else: fp.write(key+' '+key+'\n')
	#---extra tool paths
	for key in tool_paths.keys():
		fp.write(key+' '+tool_paths[key]+'\n')
	fp.close()
	
	#---write cluster-header if possible
	scratch_suffix = ''
	header_source_footer = None
	if proc_settings != None:
		if 'scratch' in proc_settings and proc_settings['scratch']: scratch_suffix = '_scratch'
	#---script repo was used to consolidate cluster_header globals in one place
	script_repo = dict([(k,globals()[k]) for k in [j for j in globals() if re.match('^cluster_header',j)]])
	if hostname != None and 'cluster_header_'+system_id+scratch_suffix in script_repo:
		scripttext = script_repo['cluster_header_'+system_id+scratch_suffix]
		#---if the header is a list, then it must contain a header and a footer
		if type(scripttext) == list:
			header_source = scripttext[0]
			header_source_footer = scripttext[1]
		else: header_source = scripttext
		header_source_mod = []
		for line in header_source.split('\n'):
			if line[:7] == '#PBS -l' and proc_settings != None:
				line = \
					'#PBS -l nodes='+str(proc_settings['nodes'])+\
					':ppn='+str(proc_settings['ppn'])+\
					(',pmem='+proc_settings['pmem'] if 'pmem' in proc_settings.keys() else '')+\
					',walltime='+str(proc_settings['walltime'])+':00:00'
			header_source_mod.append(line)
		header_source_mod = '\n'.join(header_source_mod)
	else: header_source_mod = None
	
	return proc_settings,header_source_mod,header_source_footer
	
def ultrasweep(hypothesis_default,sweep):

	"""
	Code for sweeping an arbitrarily deep dictionary over many dimensions in combinations.
	"""

	#---extract a list of lists ('t') of parameters to sweep over in combinations
	t = [i['values'] for i in sweep]
	for si,sweepset in enumerate(t):
		if type(sweepset) == dict:
			if len(sweepset.keys())==1 and 'file_list_file' in sweepset.keys():
				with open(sweepset['file_list_file'],'r') as fp:
					files = [i.strip('\n') for i in fp.readlines() if
						os.path.isfile(os.path.abspath(i.strip('\n')))]
				t[si] = files
			else: raise Exception('except: unclear sweep substitution '+str(sweepset))

	#---note that this non-numpythonic way of doing this has not been rigorously tested
	#---note that the previous meshgrid method did not work on all types
	#---we take all combinations by starting with the first list and adding all others to it
	allcombos = list([[i] for i in t[0]])
	for s in t[1:]:
		for bi in range(len(allcombos)):
			b = allcombos.pop(0)
			for r in list(s): allcombos.append(b + [r])

	#---assemble a list of hypotheses from all possible combinations of the sweep values
	#---note that this code is general, and works for an arbitrarily deep dictionary
	hypotheses = []
	#---for each combo generate a new hypothesis
	for combo in allcombos:
		#---start with the default hypothesis
		newhypo = deepcopy(hypothesis_default)
		#---each combo has a value and a route which is a sequence of dictionary keys
		#---...we loop over each route to set each final value for the sweep
		for routenum in range(len(sweep)):
			#---to get to the deepest part of that route we use tmp as a pointer
			#---...and iteratively traverse one level until the second to last level
			tmp = newhypo[sweep[routenum]['route'][0]]
			#---the following checks if we are already at the end of the dictionary 
			if type(newhypo[sweep[routenum]['route'][0]]) != dict:
				newhypo[sweep[routenum]['route'][0]] = combo[routenum]
			else:
				for i in sweep[routenum]['route'][1:-1]: tmp = tmp[i]
				#---at the final level, we now have a pointer to the lowest dictionary to set the value
				tmp[sweep[routenum]['route'][-1]] = combo[routenum]
		#---once we set all the values, the hypothesis is ready
		hypotheses.append(newhypo)	
	return hypotheses
	
def bp():

	"""
	Shortcut for a crude breakpoint.
	"""
	
	raw_input('\nbreakpoint!\n')
	
#---deep dive into a dictionary of dictionaries ad infinitum
def delve(o,*k): return delve(o[k[0]],*k[1:]) if len(k)>1 else o[k[0]]
#---check if there are any redundant elements in a list
def redundant(x): return len(set(x))!=len(x)

def mapdict(tree,routes=()):

	"""
	Converts a nested dictionary into a list of all sequences of keys which lead to terminal dictionaries.
	"""

	if (type(tree)==dict and 
		all([type(tree[i])!=dict for i in tree])):	
		yield routes
	else:
		for key in tree: 
			for route in mapdict(tree[key],routes=routes+(key,)):
				yield route

