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

helpstring = """\n
    Automacs (AMX) CONTROLLER.\n
    make script <simulation> : prepare a simulation and associated scripts 
                               see simulation types below for more details
    make script restart      : generates a CPT or GRO restart from an 
                               input-md-in.mdp in repo folder
    make upload              : prepare an upload list and rsync command for 
                               a continuation if the last step was "sim"
    make upload step=<dir>   : prepare an upload list for any step 
                               including the whole directory
    make rescript            : generate continue scripts for a simulation 
                               moved to a cluster  
    make clean               : reset the folder by deleting simulation data
                               note that the repo folder is exempt 

    simulations:
    -----------
    cgmd-bilayer         : coarse-grained bilayer under MARTINI force field
    cgmd-protein         : coarse-grained protein under MARTINI force field 
                           with structure options
    aamd-bilayer         : bilayer under CHARMM36 force field
    aamd-protein         : protein in a water box using CHARMM27 force field
    protein-homology     : single- or muliple-template homology models with 
                           options for batches of mutations
    multiply nx=N ny=M   : replicate a simulation box
    restart              : generic restart which autodetects a CPT or GRO file
    """

document_string = \
	"""
	Documentation is provided via python-Sphinx/sphinx-apidoc but is not 
	generally updated with each git push so you may	want to regenerate it.\n
	"""+str('file://'+os.path.abspath(os.path.curdir)+'/docs/_build/html/index.html')+"\n\n"

#---BASH SCRIPT CONSTANTS
#-------------------------------------------------------------------------------------------------------------

#---bash commands referenced in script_dict are recorded here
#---note that "preps" lists contain bash commands to copy general steps and connect the steps
#---note that "scripts" lists contain bash commands to run the steps in sequence

prepare_equil_sim = [
"""
#---copy equilibration files
if ! [ -d "$step_equilibration" ]; then cp -r ./sources/general-equil ./$step_equilibration
else cp -r ./sources/general-equil/* ./$step_equilibration;fi\n
#---connect equilibration settings to the previous step
sed -i 's/^SOURCEDIR.*/SOURCEDIR=\.\.\/'$step_build'/g' $step_equilibration/settings.sh
sed -i 's/cp $AMXPATH\/.*/cp $AMXPATH\/'$input_files'\/\* \.\//g' $step_equilibration/settings.sh
""",
"""
#---copy simulation continuation files
if ! [ -d "$step_simulation" ]; then cp -r ./sources/general-sim ./$step_simulation
else cp -r ./sources/general-sim/* ./$step_simulation;fi\n
#---connect simulation settings to the previous step
sed -i 's/^SOURCEDIR.*/SOURCEDIR=\.\.\/'$step_equilibration'/g' $step_simulation/settings.sh
sed -i 's/cp $AMXPATH\/.*/cp $AMXPATH\/'$input_files'\/\* \.\//g' $step_simulation/settings.sh
""",
]

prepare_restart = [
"""
#---copy simulation continuation files
mkdir ./$step_simulation
cp -r ./sources/general-sim-restart/* ./$step_simulation;
cp -r ./repo/* ./$step_simulation;
#---if detect_previous_step is set in the header then copy checkpoint files
#---...note that we use the prep scripts in a self-contained way, to avoid excessive python code
if ! [ -z "${detect_previous_step+xxx}" ]; then 
#---locate the latest checkpoint from the equilibrate folder
THISDIR=$(pwd)
cd $detect_previous_step
PRUN=0
for file in md.part*.cpt
do
	if [ $(echo ${file:7:4} | sed 's/^0*//') -gt $PRUN ]; then 
		PRUN=$(echo ${file:7:4} | sed 's/^0*//')
	fi
done
cd $THISDIR
#---copy checkpoint and input files
cp $detect_previous_step/$(printf md.part%04d.tpr $PRUN) ./$step_simulation/
cp $detect_previous_step/$(printf md.part%04d.cpt $PRUN) ./$step_simulation/
fi
"""
]

prepare_rescript = [
"""
sed -i 's/^STOPHOURS.*/STOPHOURS=\"'$walltime'\"/g' $detect_previous_step/settings.sh
"""
]

prepare_multiply = [
"""
#---copy simulation continuation files
mkdir ./$step_simulation
cp -r ./sources/general-sim-restart/* ./$step_simulation;
"""
]

call_equil = """#---execute equilibration step
cd $step_equilibration
./script-equilibrate.pl
if [[ $(tail log-script-master) =~ fail$ ]];then exit;fi
cd ..
\n"""

call_sim = """#---execute simulation (continuation) step
cd $step_simulation
./script-sim.pl
if [[ $(tail log-script-master) =~ fail$ ]];then exit;fi
cd ..
\n"""

call_sim_restart = """#---execute simulation restart step
cd $step_simulation
./script-sim-restart.pl
if [[ $(tail log-script-master) =~ fail$ ]];then exit;fi
cd ..
\n"""

call_sim_rescript = """#---execute simulation restart step
cd $detect_previous_step
./script-sim.pl
if [[ $(tail log-script-master) =~ fail$ ]];then exit;fi
cd ..
\n"""

call_sim_multiply = """#---execute simulation restart step
cd $step_simulation
./script-sim-restart.pl
if [[ $(tail log-script-master) =~ fail$ ]];then exit;fi
cd ..
\n"""

python_from_bash = """#---execute python step
python -c "execfile('./amx/header');
try: COMMAND;
except Exception,e:
	print str(e);
	traceback.print_exc();
	print 'fail';"
go=$(tail ROOTDIR/log-script-master)
if [[ "$go" =~ fail$ ]];then exit;fi
\n"""

#---BASH SCRIPT GENERATOR
#-------------------------------------------------------------------------------------------------------------

#---amx procedure details
script_dict = {
	'cgmd-bilayer':{
		'prep':[prepare_equil_sim[0]+\
			"sed -i 's/STEPSTRING.*/STEPSTRING=\"nvt npt\"/g' $step_equilibration/settings.sh"]+\
			prepare_equil_sim[1:],
		'steps':{
			'step_monolayer':'s01-build-lipidgrid',
			'step_build':'s02-build-bilayer',
			'step_equilibration':'s03-equil',
			'step_simulation':'s04-sim',
			'input_files':'cgmd-bilayer-equil',
			},
		'continue':True,
		'sequence':[
			multiresub(dict({
				'ROOTDIR':'s01-build-lipidgrid',
				'COMMAND':'MonolayerGrids(rootdir=\'s01-build-lipidgrid\')',
				}),python_from_bash),
			multiresub(dict({
				'ROOTDIR':'s02-build-bilayer',
				'COMMAND':'Bilayer(rootdir=\'s02-build-bilayer\',previous_dir=\'s01-build-lipidgrid\')',
				}),python_from_bash),
			call_equil,
			call_sim,
			],
		},
	'aamd-bilayer':{
		'prep':prepare_equil_sim,
		'steps':{
			'step_monolayer':'s1-build-lipidgrid',
			'step_build':'s02-build-bilayer',
			'step_equilibration':'s03-equil',
			'step_simulation':'s04-sim',
			'input_files':'aamd-bilayer-equil',
			},
		'continue':True,
		'sequence':[
			multiresub(dict({
				'ROOTDIR':'s01-build-lipidgrid',
				'COMMAND':'MonolayerGrids(rootdir=\'s01-build-lipidgrid\')',
				}),python_from_bash),
			multiresub(dict({
				'ROOTDIR':'s02-build-bilayer',
				'COMMAND':'Bilayer(rootdir=\'s02-build-bilayer\',previous_dir=\'s01-build-lipidgrid\')',
				}),python_from_bash),
			call_equil,
			call_sim,
			],
		},
	'aamd-protein':{
		'prep':[prepare_equil_sim[0]+\
			"sed -i 's/STEPSTRING.*/STEPSTRING=\"nvt npt\"/g' $step_equilibration/settings.sh"]+\
			prepare_equil_sim[1:],
		'steps':{
			'step_build':'s01-build-protein-water',
			'step_equilibration':'s02-equil',
			'step_simulation':'s03-sim',
			'input_files':'aamd-protein-equil',
			'detect_previous_step':None,
			},
		'continue':True,
		'sequence':[
			multiresub(dict({
				'ROOTDIR':'$step_build',
				'COMMAND':'ProteinWater(rootdir=\'$step_build\')',
				}),python_from_bash),
			call_equil,
			call_sim,		
			],
		},
	'cgmd-protein':{
		'prep':[prepare_equil_sim[0]+\
			"sed -i 's/STEPSTRING.*/STEPSTRING=\"nvt-short nvt npt\"/g' $step_equilibration/settings.sh"]+\
			prepare_equil_sim[1:],
		'steps':{
			'step_build':'s01-build-protein-water',
			'step_equilibration':'s02-equil',
			'step_simulation':'s03-sim',
			'input_files':'cgmd-protein-equil',
			},
		'continue':True,
		'sequence':[
			multiresub(dict({
				'ROOTDIR':'s01-build-protein-water',
				'COMMAND':'ProteinWater(rootdir=\'s01-build-protein-water\')',
				}),python_from_bash),
			call_equil,
			call_sim,		
			],
		},
	'protein-homology':{
		'prep':[],
		'steps':{
			'step_homology':'s01-homology',
			'input_files':'aamd-protein-homology',
			},
		'sequence':[
			multiresub(dict({
				'ROOTDIR':'s01-homology',
				'COMMAND':'ProteinHomology(rootdir=\'s01-homology\')',
				}),python_from_bash),
			],
		},
	'restart':{
		'prep':prepare_restart,
		'steps':{
			'step_simulation':'s01-restart',
			'detect_previous_step':None,
			},
		'continue':True,
		'sequence':[
			call_sim_restart,
			],
		},
	'rescript':{
		'prep':prepare_rescript,
		'steps':{
			'detect_previous_step':None,
			},
		'sequence':[
			call_sim_rescript,
			],
		},
	'multiply':{
		'prep':prepare_multiply,
		'steps':{
			'step_simulation':'s01-multiply',
			'detect_previous_step':None,
			},
		'continue':True,
		'sequence':[
			multiresub(dict({
				'ROOTDIR':'$step_simulation',
				'COMMAND':'Multiply(rootdir=\'$step_simulation\',previous_dir=\'$previus_step\','+\
					'nx=\'$nx\',ny=\'$ny\',nz=\'$nz\')',
				}),python_from_bash),
			call_sim_restart,
			call_sim,
			],
		},
	'cgmd-protein-bilayer':{
		'prep':[prepare_equil_sim[0]+\
			"sed -i 's/STEPSTRING.*/STEPSTRING=\"npt\"/g' $step_equilibration/settings.sh"]+\
			prepare_equil_sim[1:],
		'steps':{
			'step_build':'s01-protein-bilayer',
			'input_files':'cgmd-protein-bilayer',
			'step_equilibration':'s02-equil',
			'step_simulation':'s03-sim',
			},
		'sequence':[
			multiresub(dict({
				'ROOTDIR':'$step_build',
				'COMMAND':'ProteinBilayer(rootdir=\'$step_build\')',
				}),python_from_bash),
			call_equil,
			call_sim,
			],
		},
	}

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def clean(sure=True,protected=False,extras=None):
	'''
	Clears all run files from the current directory structure.
	'''
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
					
def docs(extras=None):
	'''Regenerate the documentation using sphinx-apidoc and code in amx/gendoc.sh'''
	os.system('./amx/script-make-docs.sh '+os.path.abspath('.'))
	
def upload(step=None,extras=None,silent=False):
	'''Provides a file list and rsync syntax to upload this to a cluster without extra baggage.'''
	gofiles = ['controller.py','makefile','settings']
	startstep,oldsteps = chain_steps()
	if step != None:
		if step not in oldsteps: raise Exception('except: cannot find that step folder')
		last = step
	else: last = oldsteps[-1]
	if not silent: print last
	for root,dirnames,filenames in os.walk('./amx'): break
	for fn in filenames: gofiles.append('amx/'+fn)
	if step == None:
		for t in latestcheck(last): gofiles.append(last+'/'+t)
	else:
		for root,dirnames,filenames in os.walk(last): break
		for fn in filenames: gofiles.append(last+'/'+fn)
	for root,dirnames,filenames in os.walk(last): break
	for j in [i for i in filenames if re.match('.+\.(pl|sh)$',i)]: gofiles.append(last+'/'+j)
	cwd = os.path.basename(os.path.abspath(os.getcwd()))
	with open('upload-rsync-list.txt','w') as fp: 
		for i in gofiles: fp.write(i+'\n')
	if not silent: print 'upload list written to upload-rsync-list.txt\nrun '+\
		'the following rsync to your destination'
	if not silent: print 'rsync -avi --files-from=upload-rsync-list.txt ../'+cwd+' DEST/'+cwd
	
def copycode():
	startstep,oldsteps = chain_steps()
	cwd = os.path.basename(os.path.abspath(os.getcwd()))
	with open('upload-rsync-code-excludes.txt','w') as fp: 
		for step in oldsteps: fp.write(step+'\n')
	print 'rsync -avi --exclude-from=upload-rsync-code-excludes.txt ./ ../automacs-code'

#-------------------------------------------------------------------------------------------------------------

def batchscript(batchdir,proc):
	'''
	Recurseively run make commands for a large batch.
	'''
	batchdir = os.path.abspath(batchdir)+'/'
	for root,dirnames,filenames in os.walk(batchdir,followlinks=False): break
	for dn in dirnames: 
		print batchdir+dn
		print 'calling = '+'python controller.py script '+proc
		call('python controller.py script '+proc,cwd=batchdir+dn,suppress_stdout=True)

def rescript(single=None,**extras): script(single=single,rescript=True,**extras)
					
def script(single=None,rescript=False,**extras):

	'''
	Deposit PBS scripts in each step folder for the corresponding step and make master scripts.
	'''
	
	#---check for sensible script target
	if single == None and rescript == False: raise Exception('except: unclear target')
	#---a rescript on the cluster will find the most recent step and add a continue script to it
	elif single == None and rescript == True:
		startstep,oldsteps = chain_steps()
		proc_settings,header_source_mod,hs_footer = get_proc_settings()
		if header_source_mod == None: raise Exception('except: lone rescript only works on clusters')
		#---run prep scripts
		extras['walltime'] = proc_settings['walltime']
		prep_scripts('rescript',script_dict,extras=extras)
		#---write a list of files to move to scratch
		upload(silent=True)
		with open('upload-rsync-list.txt','a') as fp: fp.write('gmxpaths.conf\n')
		#---write a simulation continuation script in the root directory on clusters
		fp = open('cluster-md-continue','w')
		fp.write(header_source_mod)
		fp.write(script_maker('rescript',script_dict,
			module_commands=proc_settings['module'],
			sim_only=None,extras=extras,))
		if hs_footer != None: fp.write(hs_footer)
		fp.close()		
		
	#---single procedure script maker
	else:
		target = single
		print helpstring
		print "\tGenerating singleton script.\n"

		#---update the input specs file simscale parameter according to the target simulation
		#---? note that this might be broken
		if any([j=='cgmd' for j in target.split('-')]): simscale = 'cgmd'                                                                                                       
		elif any([j=='aamd' for j in target.split('-')]): simscale = 'aamd'                                                                                                     
		else: simscale = None                                                                                                                                               
		if simscale != None:                                                                                                                                                    
			for fn in glob.glob('./inputs/input-specs*'):                                                                                                                   
				sub = subprocess.call(['sed','-i',"s/^.*simscale.*$/simscale = \'"+simscale+"\'/",fn])
	
		proc_settings,header_source_mod,hs_footer = get_proc_settings()

		#---write the local script
		print '\twriting script: script-master-'+str(target)
		fp = open('script-master-'+str(target),'w')
		fp.write('#!/bin/bash\n')
		fp.write(script_maker(target,script_dict,extras=extras))
		fp.close()
		os.system('chmod u+x '+'script-master-'+str(target))
		
		#---prepare the files and directories
		if not rescript: prep_scripts(target,script_dict,extras=extras)
		
		#---for simulations that can be continued we write a script for only the final step in the sequence
		if 'continue' in script_dict[target].keys() and script_dict[target]['continue']:
			#---use chain_steps to find the final folder to deposit a continuation script if necessary
			startstep,oldsteps = chain_steps()
			fp = open(oldsteps[-1]+'/script-md-continue','w')
			fp.write('#!/bin/bash\n')
			fp.write(script_maker(target,script_dict,sim_only=call_sim,extras=extras))
			fp.close()
			os.system('chmod u+x '+oldsteps[-1]+'/script-md-continue')

		print '\texecute locally with ./'+'script-master-'+str(target)
		if rescript: print '\tsince this is a rescript you probably want to use the continue script'
		print '\tsee the documentation for details\n'
		
def batch(**extras):
	'''
	Function which spawns a particular simulation into many.\n
	Recommended useage: 
		make clean sure; make script protein-homology; ./script-master-protein-homology
		make spawn proc=aamd-protein infiles=repo/simulation_targets.txt batchdir=../protein-v1020-batch
	'''
	batchspecs = {}
	execfile('inputs/input-specs-batch.dat',batchspecs)
	batchdir = os.path.abspath(batchspecs['batchdir'])+'/'
	if os.path.isdir(batchdir): raise Exception('except: batch directory already exists')
	else: os.mkdir(batchdir)
		
	#---get the default specs file
	specsfile = {
		'aamd-protein':'inputs/input-specs-protein.dat',
		'cgmd-protein':'inputs/input-specs-protein.dat',
		'aamd-bilayer':'inputs/input-specs-bilayer.dat',
		'cgmd-protein':'inputs/input-specs-bilayer.dat',
		}[batchspecs['procedure']]
	defaults = {}
	execfile(specsfile,defaults)
	del defaults['__builtins__']
	
	#---set any overrides in the default dictionary
	overrides = batchspecs['overrides']
	for key in overrides:
		if key in defaults.keys(): defaults[key] = overrides[key]
		else: raise Exception('except: overrides key not found in default input-specs')
	
	#---create the sweep
	hypotheses = ultrasweep(defaults,batchspecs['sweep'])
	
	#---for each simulation copy automacs, input_filename, and input-specs
	loclist = []
	for hi,hypo in enumerate(hypotheses):
		newdir = batchdir+'sim-v'+batchspecs['callsign']+'-v'+('%03d'%(hi+1))+'/'
		os.mkdir(newdir)
		os.mkdir(newdir+'/inputs')
		with open(newdir+specsfile,'w') as fp:
			fp.write('#!/usr/bin/python\n#---autogenerated specs file\n')
			for key in hypo:
				keystring = json.dumps(hypo[key],indent=4).replace('false','False').replace('true','True')
				fp.write(key+' = '+keystring+'\n')
		if 'input_filename' in hypo.keys():
			#---recursively create directories to copy the full path of the input_filename
			dirlist = hypo['input_filename'].strip('./').split('/')
			for di,dn in enumerate(dirlist):
				if di>0: os.mkdir(newdir+'/'.join(dirlist[:di]))
			copy(hypo['input_filename'],newdir+hypo['input_filename'])
			#---repo is excluded for size concerns but can be added if necessary
			for fn in ['amx','controller.py','sources','settings','makefile','repo'][:-1]:
				copy('./'+fn,newdir+fn)
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
	

#---INTERFACE
#-------------------------------------------------------------------------------------------------------------

'''
Functions exposed to makefile:
def avail
def timeslice
def catalog
'''

def makeface(arglist):

	'''
	Interface to makefile.
	'''
	
	#---print help if no arguments
	if arglist == []:
		print niceblock(helpstring,newlines=True)
		return
		
	#---we prepare a kwargs and an args variable to send to the next function
	#---note that we generally just use kwargs exclusively for clarity
	kwargs = dict()
	
	#---always get the function name from the first argument
	func = arglist.pop(0)
	#---define the arguments expected for each function
	#---note that 'functionname':{'args':[],'module_name':None} is the default empty interface to a python
	#---...function however you must remember to add keywords to args if you expect them to be standalone
	#---...if you add keywords to args they can be True (if passed on the command line) otherwise defaulted
	#---...to false while any var=val commands will be passed along via extras
	argdict = {
		'clean':{'args':['protected','sure'],'module_name':None,'defaults':{'sure':False}},
		'script':{'module_name':None,'defaults':[],'args':[],
			'singles':[
				'aamd-bilayer',
				'cgmd-bilayer',
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
		if not re.match('^[a-z,A-Z,0-9]+\=([a-z,A-Z,0-9,\-,_,\.,\/]+)+$',a):
			kwargs[a] = True if a in argd['args'] else False
			if a in argd['args']: arglist.remove(a)
		else:
			flag = re.findall('^([a-z,A-Z,0-9]+)\=([a-z,A-Z,0-9,\-,_,\.,\/]+)$',a)[0]
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
	
