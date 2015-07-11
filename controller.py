#!/usr/bin/python

import os,re,sys
import subprocess
import socket
from copy import deepcopy
import glob
import json
import numpy as np
from amx.tools import copy
from amx.tools import call

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

execfile('settings')
from amx.tools import checkout,multiresub,confirm
from amx.tools import script_maker,prep_scripts,niceblock,argsort,chain_steps,latestcheck
from amx.tools import get_proc_settings,ultrasweep

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def clean(sure=True,protected=False,extras=None):

	"""
	Clears all run files from the current directory structure.
	"""

	if protected: print 'delete protected locations = '+str(protected)
	rootdir = os.path.expanduser('.')
	delfiles,delfiles_protect,deldirs,deldirs_protect,deldirs_all = [],[],[],[],[]
	for root,dirnames,filenames in os.walk(rootdir,followlinks=True):
		if protected:
			if re.match('^\.?\/?(sources[^\/]|repo)',root):
				deldirs_protect.append([root,''])
		if re.match('^\.?\/?[a-z][0-9]{1,2}-',root):
			deldirs.append([root,''])
		for filename in filenames:
			if any([re.match(i,filename) for i in patlist]) and \
				not any([re.match(w,root) for w in whitelist]) and \
				not any([re.match(w,filename) for w in whitelist]) and \
				not re.search('source',root) and \
				not re.search('input',root) and \
				not re.search('repo',root) and \
				not re.search('amx',root):
				delfiles.append([root,filename])
			if protected:
				if re.match('^\.?\/?(sources[^\/]|repo).*\/?.+',root+filename):
					delfiles_protect.append([root,filename])
	if any([i != [] for i in [delfiles,delfiles_protect,deldirs,deldirs_protect,deldirs_all]]):
		if delfiles != []: print 'files to delete:'
		for i in delfiles:  print ' '+i[0]+'/'+i[1]
		if delfiles_protect != []: print 'protected files to delete:'
		for i in delfiles_protect: print ' '+i[0]+'/'+i[1]
		if deldirs != []: print 'directories to delete:'
		for i in deldirs:  print ' '+i[0]+'/'+i[1]
		if deldirs_protect != []: print 'protected directories to delete:'
		for i in deldirs_protect: print ' '+i[0]+'/'+i[1]
		if deldirs_all != []: print 'complete directories to delete:'
		for i in deldirs_all: print ' '+i[0]+'/'+i[1]
		if not sure and not confirm(): 
			print 'cancel'
			return
		for f in delfiles: os.system('rm '+f[0]+'/'+f[1])
		for f in delfiles_protect: os.system('rm '+f[0]+'/'+f[1])
		for f in deldirs:
			if [re.match(i,f[0]) for i in folder_delete]: os.system('rm -rf '+f[0]+'/'+f[1])
			else: os.system('rmdir '+f[0]+'/'+f[1])
		for f in deldirs_protect: os.system('rm -rf '+f[0]+'/'+f[1])
					
def docs(clean=False):

	"""
	Regenerate the documentation using sphinx-apidoc and code in amx/gendoc.sh
	"""

	import shutil
	if clean: shutil.rmtree('./docs')
	else:
		if not os.path.isfile('gmxpaths.conf'): 
			with open('gmxpaths.conf','w') as fp: fp.write('# placeholder')
		os.system('./amx/script-make-docs.sh '+os.path.abspath('.'))
	
def upload(step=None,extras=None,silent=False,scriptfile=None):

	"""
	Provides a file list and rsync syntax to upload this to a cluster without extra baggage.
	"""

	gofiles = ['controller.py','makefile','settings']
	startstep,oldsteps = chain_steps()
	if step != None:
		if step not in oldsteps: raise Exception('except: cannot find that step folder')
		last = step
	else: last = oldsteps[-1]
	for root,dirnames,filenames in os.walk('./amx'): break
	for fn in filenames: gofiles.append('amx/'+fn)
	if step == None:
		lcheck = latestcheck(last)
		if lcheck == []: return False
		for t in latestcheck(last): gofiles.append(last+'/'+t)
	else:
		for root,dirnames,filenames in os.walk(last): break
		for fn in filenames: gofiles.append(last+'/'+fn)
	for root,dirnames,filenames in os.walk(last): break
	for j in [i for i in filenames if re.match('.+\.(pl|sh)$',i)]: gofiles.append(last+'/'+j)
	cwd = os.path.basename(os.path.abspath(os.getcwd()))
	with open('upload-rsync-list.txt','w') as fp: 
		for i in gofiles: fp.write(i+'\n')
	#---custom functionality only used by batch (so be careful with relative paths)
	if scriptfile != None:
		print "TRYING TO OPEN "+scriptfile
		with open(scriptfile,'a') as fp:
			fp.write('rsync -avi --files-from='+cwd+'/upload-rsync-list.txt '+cwd+' DEST/'+cwd+'\n')
	else:
		if not silent: print 'upload list written to upload-rsync-list.txt\nrun '+\
			'the following rsync to your destination'
		if not silent: print 'rsync -avi --files-from=upload-rsync-list.txt ../'+cwd+' DEST/'+cwd
	
def copycode():

	"""
	Utility function for moving code. Deprecated.
	"""

	startstep,oldsteps = chain_steps()
	cwd = os.path.basename(os.path.abspath(os.getcwd()))
	with open('upload-rsync-code-excludes.txt','w') as fp: 
		for step in oldsteps: fp.write(step+'\n')
	print 'rsync -avi --exclude-from=upload-rsync-code-excludes.txt ./ ../automacs-code'

#-------------------------------------------------------------------------------------------------------------

def batchscript(batchdir,proc):

	"""
	Recurseively run make commands for a large batch.
	"""

	batchdir = os.path.abspath(batchdir)+'/'
	for root,dirnames,filenames in os.walk(batchdir,followlinks=False): break
	for dn in dirnames: 
		print batchdir+dn
		print 'calling = '+'python controller.py script '+proc
		call('python controller.py script '+proc,cwd=batchdir+dn,suppress_stdout=True)

def rescript(single=None,**extras): script(single=single,rescript=True,**extras)
					
def script(single=None,rescript=False,**extras):

	"""
	Deposit PBS scripts in each step folder for the corresponding step and make master scripts.
	"""
	
	#---flags which are passed to settings.sh files if necessary
	sets_flags = ['gpu_flag']
	sets_pass = dict([(sf,extras[sf]) for sf in [s for s in extras if s in sets_flags]])

	#---check for sensible script target
	if single == None and rescript == False: raise Exception('except: unclear target')
	#---a rescript on the cluster will find the most recent step and add a continue script to it
	elif single == None and rescript == True:
		startstep,oldsteps = chain_steps()
		proc_settings,header_source_mod,hs_footer = get_proc_settings()
		if header_source_mod == None: raise Exception('except: lone rescript only works on clusters')
		#---run prep scripts
		#---note that we add a 2% buffer here combined with 1% on the maxh flag so enough time to copy
		#---...files back from scratch
		extras['walltime'] = proc_settings['walltime']*0.98
		prep_scripts('rescript',script_dict,extras=extras,extra_settings=sets_pass)
		#---write a list of files to move to scratch
		#---use the upload routine to only copy necessary files to scratch
		if not ('start' in extras.keys() and extras['start']):
			upload(silent=True)
			with open('upload-rsync-list.txt','a') as fp: fp.write('gmxpaths.conf\n')
		#---if start flag is present then move everything
		else:
			files = [os.path.join(dp, f) for dp, dn, fns in os.walk('./') for f in fns]
			files = [re.findall('^'+'\.?\/?'+'((?!\.git).+)',f) for f in files]
			files = [f[0] for f in files if f!=[]]
			with open('upload-rsync-list.txt','w') as fp:
				for f in files: fp.write(f+'\n')

		#---write a simulation continuation script in the root directory on clusters
		fp = open('cluster-md-continue','w')
		fp.write(header_source_mod)
		fp.write(script_maker('rescript',script_dict,
			module_commands=proc_settings['module'],
			sim_only=None,extras=extras))
		if hs_footer != None: fp.write(hs_footer)
		fp.close()		
		
	#---single procedure script maker
	#---note that supplying rescript=True will skip the preparation scripts and remake cluster scripts
	else:
		target = single
		print helpstring
		print "\tGenerating singleton script.\n"
		if rescript: 
			if extras == None: extras = {'rescript':True}
			else: extras['rescript'] = True

		#---update the input specs file simscale parameter according to the target simulation
		#---? note that this might be broken
		if any([j=='cgmd' for j in target.split('-')]): simscale = 'cgmd'
		elif any([j=='aamd' for j in target.split('-')]): simscale = 'aamd'
		else: simscale = None
		if simscale != None:
			for fn in glob.glob('./inputs/input_specs*'):
				sub = subprocess.call(['sed','-i',"s/^.*simscale.*$/simscale = \'"+simscale+"\'/",fn])
		proc_settings,header_source_mod,hs_footer = get_proc_settings()
		#---if possible write the cluster script
		if header_source_mod != None or rescript:
			with open('./cluster-'+target,'w') as fp:
				fp.write(header_source_mod)
				fp.write(script_maker(target,script_dict,module_commands=proc_settings['module'],
					sim_only=None,extras=extras))
				if hs_footer != None: fp.write(hs_footer)
				if 'carefultime' in extras.keys() and extras['carefultime']: 
					fp.write('qsub cluster-md-continue')
				#---use the upload routine to only copy necessary files to scratch
				extras_start = not any([re.match('^[a-z][0-9]{1,2}',i) for i in os.listdir('.')])
				if not ('start' in extras.keys() and extras['start']) and not extras_start:
					upload(silent=True)
					with open('upload-rsync-list.txt','w') as fp: fp.write('gmxpaths.conf\n')
				#---if start flag is present then move everything
				else:
					files = [os.path.join(dp, f) for dp, dn, fns in os.walk('./') for f in fns]
					files = [re.findall('^'+'\.?\/?'+'((?!\.git).+)',f) for f in files]
					files = [f[0] for f in files if f!=[]]
					with open('upload-rsync-list.txt','w') as fp:
						for f in files: fp.write(f+'\n')
		
		#---write the local script
		print '\twriting script: script-'+str(target)
		fp = open('script-'+str(target),'w')
		fp.write('#!/bin/bash\n')
		fp.write(script_maker(target,script_dict,extras=extras))
		fp.close()
		os.system('chmod u+x '+'script-'+str(target))
		
		#---prepare the files and directories
		if not rescript: prep_scripts(target,script_dict,extras=extras,extra_settings=sets_pass)
		
		#---for simulations that can be continued we write a script for only the final step in the sequence
		if 'continue' in script_dict[target].keys() and script_dict[target]['continue']:
			#---use chain_steps to find the final folder to deposit a continuation script if necessary
			startstep,oldsteps = chain_steps()
			fp = open(oldsteps[-1]+'/script-md-continue','w')
			fp.write('#!/bin/bash\n')
			fp.write(script_maker(target,script_dict,sim_only=call_sim,extras=extras))
			fp.close()
			os.system('chmod u+x '+oldsteps[-1]+'/script-md-continue')

		print '\texecute locally with ./'+'script-'+str(target)
		if rescript: print '\tsince this is a rescript you probably want to use the continue script'
		print '\tsee the documentation for details\n'
		
def batch(**extras):

	"""
	Function which spawns a particular simulation into many.\n
	Recommended useage: 
		make clean sure; make script protein-homology; ./script-protein-homology
		make spawn proc=aamd-protein infiles=repo/simulation_targets.txt batchdir=../protein-v1020-batch
	"""

	batchspecs = {}
	execfile('inputs/input_specs_batch.py',batchspecs)
	batchdir = os.path.abspath(batchspecs['batchdir'])+'/'
	if os.path.isdir(batchdir): raise Exception('except: batch directory already exists')
	else: os.mkdir(batchdir)
		
	#---get the default specs file
	specsfile = {
		'aamd-protein':'inputs/input_specs_protein.py',
		'cgmd-protein':'inputs/input_specs_protein.py',
		'aamd-bilayer':'inputs/input_specs_bilayer.py',
		'cgmd-protein':'inputs/input_specs_bilayer.py',
		'cgmd-bilayer':'inputs/input_specs_bilayer.py',
		}[batchspecs['procedure']]
	defaults = {}
	execfile(specsfile,defaults)
	del defaults['__builtins__']
	
	#---set any overrides in the default dictionary
	overrides = batchspecs['overrides']
	for key in overrides:
		if key in defaults.keys(): defaults[key] = overrides[key]
		else: raise Exception('except: overrides key not found in default input_specs')
	
	#---create the sweep
	hypotheses = ultrasweep(defaults,batchspecs['sweep'])
	
	#---for each simulation copy automacs, input_filename, and input_specs
	loclist = []
	for hi,hypo in enumerate(hypotheses):
		newdir = batchdir+'sim-v'+batchspecs['callsign']+'-v'+('%03d'%(hi+1))+'/'
		os.mkdir(newdir)
		os.mkdir(newdir+'/inputs')
		if 'input_filename' in hypo.keys():
			#---recursively create directories to copy the full path of the input_filename
			dirlist = hypo['input_filename'].strip('./').split('/')
			#---we recreate the full directory structure of the source files
			dir_start_level = dirlist.index(os.path.basename(os.getcwd()))
			for di,dn in enumerate(dirlist):
				if di>dir_start_level+1: os.mkdir(newdir+'/'.join(dirlist[dir_start_level+1:di]))
			newdirpath = '/'.join(dirlist[dir_start_level+1:])
			hypo['input_filename'] = newdirpath
			copy(hypo['input_filename'],newdir+hypo['input_filename'])
			#---repo is excluded for size concerns but can be added if necessary
			for fn in ['amx','controller.py','sources','settings','makefile','repo'][:-1]:
				copy('./'+fn,newdir+fn)
		with open(newdir+specsfile,'w') as fp:
			fp.write('#!/usr/bin/python\n#---autogenerated specs file\n')
			for key in hypo:
				keystring = json.dumps(hypo[key],indent=4).\
					replace('false','False').\
					replace('true','True').\
					replace('null','None')
				fp.write(key+' = '+keystring+'\n')
		copy('./inputs/input_specs_mdp.py',newdir+'inputs/input_specs_mdp.py')
		call('python controller.py script aamd-protein',cwd=os.path.abspath(newdir),
			silent=True,suppress_stdout=True)
		loclist.append(newdir)
	
	#---write a summary dictionary
	hypodict = {}
	for hi,hypo in enumerate(hypotheses): hypodict[loclist[hi]] = deepcopy(hypotheses[hi])
	with open(batchspecs['batchdir']+'/batchdict.py','w') as fp:
		fp.write('#!/usr/bin/python\n#---batch simulation table of contents\n')
		keystring = json.dumps(hypodict,indent=4).replace('false','False').replace('true','True')
		fp.write('batch_toc = '+keystring+'\n')
	print 'Completed batch generation in '+batchspecs['batchdir']
	print 'Wrote a dictionary of all simulations to '+batchspecs['batchdir']+'/batchdict.py'
	with open(batchspecs['batchdir']+'script-batch-serial.sh','w') as fp:
		fp.write('#!/bin/bash\n#---batch serial execution\n')
		for key in hypodict:
			fp.write('cd '+key+'\n')
			fp.write('./script-aamd-protein\n')
			fp.write('cd ..'+'\n')
	print 'Wrote a serial execution script to script-batch-serial.sh'
	with open(batchspecs['batchdir']+'script-batch-upload.sh','w') as fp:
		fp.write('#!/bin/bash\n#---batch uploads\n')
	for key in hypodict:
		call('make upload step=s04-sim scriptfile='+
			os.path.abspath(batchspecs['batchdir']+'script-batch-upload.sh'),cwd=key)
	print 'Wrote a serial upload script to script-batch-serial.sh\n'+\
		'To use the uploader, you have to find-replace DEST with your target system\n'+\
		'from '+batchspecs['batchdir']+' run:\n'+\
		'sed -i \'s/DEST/compbio:~/g\' script-batch-upload.sh'+\
		'Wait until step four is complete before uploading.'

#---INTERFACE
#-------------------------------------------------------------------------------------------------------------

"""
Functions exposed to makefile:
def avail
def timeslice
def catalog
"""

def makeface(arglist):

	"""
	Interface to makefile.
	"""
	
	#---print help if no arguments
	if arglist == []:
		print niceblock(helpstring,newlines=True)
		return
		
	#---we prepare a kwargs and an args variable to send to the next function
	#---note that we generally just use kwargs exclusively for clarity
	kwargs = dict()
	
	#---always get the function name from the first argument
	func = arglist.pop(0)
	for stray in ['--','w','i']:
		if stray in arglist: arglist.remove(stray)
	
	#---define the arguments expected for each function
	#---note that 'functionname':{'args':[],'module_name':None} is the default empty interface to a python
	#---...function however you must remember to add keywords to args if you expect them to be standalone
	#---...if you add keywords to args they can be True (if passed on the command line) otherwise defaulted
	#---...to false while any var=val commands will be passed along via extras
	argdict = {
		'clean':{'args':['protected','sure'],'module_name':None,'defaults':{'sure':False}},
		'script':{'module_name':None,'defaults':[],'args':['carefultime','rescript','start','gpu_flag'],
			'singles':[
				'aamd-bilayer',
				'cgmd-bilayer',
				'cgmd-bilayer-sculpt',
				'aamd-protein',
				'cgmd-protein',
				'cgmd-protein-bilayer',
				'protein-homology',
				'restart',
				'multiply',
				]},
		'upload':{'args':[],'module_name':None},
		'batch':{'args':[],'module_name':None},
		'batchscript':{'args':[],'module_name':None},
		'docs':{'args':['clean'],'module_name':None},
		}
	#---rescript is an alias for script for only doing the continue scripts
	argdict['rescript'] = deepcopy(argdict['script'])
	
	#---make a copy of the dictionary for pruning
	if func == 'gitpush': return
	elif func not in argdict.keys():
		globals()[func]()
		return
	else: argd = deepcopy(argdict[func])
	
	#---special keywords are handled below
	if func == 'clean':
		for default in argd['defaults'].keys():
			kwargs[default] = argd['defaults'][default]
		for a in list(arglist): 
			kwargs[a] = True if a in argd['args'] else False
			if a in argd['args']: arglist.remove(a)
	elif func == 'script' or (func == 'rescript' and len(arglist)>0):
		if arglist[0] in argd['singles']: kwargs['single'] = arglist.pop(0)
	#---all remaining keywords are handled as flags
	for a in list(arglist):
		if not re.match('^[a-z,A-Z,0-9,_]+\=([a-z,A-Z,0-9,\-,_,\.,\/]+)+$',a):
			kwargs[a] = True if a in argd['args'] else False
			if a in argd['args']: arglist.remove(a)
		else:
			flag = re.findall('^([a-z,A-Z,0-9,_]+)\=([a-z,A-Z,0-9,\-,_,\.,\/]+)$',a)[0]
			kwargs[flag[0]] = flag[1]
			arglist.remove(a)
	if arglist != []: 
		print 'warning: unprocessed arguments: '+str(arglist)
		kwargs['extras'] = arglist

	#---call the function possibly in another module
	if argd['module_name'] == None and func != 'gitpush': target = globals()[func]
	else: target = getattr(getattr(controller,argd['module_name']),func)
	target(**kwargs)
	return	

#---MAIN
#-------------------------------------------------------------------------------------------------------------

if __name__ == "__main__": makeface(sys.argv[1:])
	
