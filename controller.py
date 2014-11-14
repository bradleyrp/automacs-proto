#!/usr/bin/python -i

import os,re,sys
import subprocess
import socket
import copy
import glob

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

execfile('settings')
from amx.tools import checkout,multiresub,confirm
from amx.tools import script_maker,prep_scripts,niceblock,argsort,chain_steps

helpstring = """\n
	Automacs (AMX) CONTROLLER.\n
	make <operation>
	tools for running automatic GROMACS simulations\n
	operations:
	----------
	clean           : clears the subdirectories of all run files
	script <method> : generates PBS and bash scripts to construct simulations
	docs            : re-generates documentation (development only)\n
	methods:
	-------
	cgmd-bilayer    : build a coarse-grained lipid bilayer
	aamd-bilayer    : build an atomistic lipid bilayer
	aamd-protein    : build a protein-in-water simulation
	"""

document_string = \
	"""
	Documentation is provided via python-Sphinx/sphinx-apidoc which
	automatically generates python-style documentation web pages 
	according to the docstrings. This means that the code is also 
	heavily commented and should be readable as well. Visit the 
	following link to view the documentation.\n
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
sed -i 's/SIMTYPE.*/SIMTYPE="'$simtype'"/g' $step_equilibration/settings.sh
sed -i 's/cp $AMXPATH\/.*/cp $AMXPATH\/'$input_files'\/\* \.\//g' $step_equilibration/settings.sh
""",
"""
#---copy simulation continuation files
if ! [ -d "$step_simulation" ]; then cp -r ./sources/general-sim ./$step_simulation
else cp -r ./sources/general-sim/* ./$step_simulation;fi\n
#---connect simulation settings to the previous step
sed -i 's/^SOURCEDIR.*/SOURCEDIR=\.\.\/'$step_equilibration'/g' $step_simulation/settings.sh
sed -i 's/SIMTYPE.*/SIMTYPE="'$simtype'"/g' $step_simulation/settings.sh
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
		'prep':prepare_equil_sim,
		'steps':{
			'step_monolayer':'s01-build-lipidgrid',
			'step_build':'s02-build-bilayer',
			'step_equilibration':'s03-equil',
			'step_simulation':'s04-sim',
			'simtype':'cgmd-bilayer',
			'input_files':'cgmd-bilayer-equil',
			},
		'continue':True,
		'sequence':['\n',
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
			'simtype':'aamd-bilayer',
			'input_files':'aamd-bilayer-equil',
			},
		'continue':True,
		'sequence':['\n',
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
			'simtype':'aamd-protein',
			'input_files':'aamd-protein-equil',
			},
		'continue':True,
		'sequence':['\n',
			multiresub(dict({
				'ROOTDIR':'s01-build-protein-water',
				'COMMAND':'ProteinWater(rootdir=\'s01-build-protein-water\')',
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
			'simtype':'cgmd-protein',
			'input_files':'cgmd-protein-equil',
			},
		'continue':True,
		'sequence':['\n',
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
		'sequence':['\n',
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
		'sequence':['\n',
			call_sim_restart,
			],
		},
	}

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def clean(sure=True,protected=False):
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
		
					
def docs():
	'''Regenerate the documentation using sphinx-apidoc and code in amx/gendoc.sh'''
	os.system('./amx/script-make-docs.sh '+os.path.abspath('.'))


#---MAKE
#-------------------------------------------------------------------------------------------------------------
					
def script(single=None):
	'''
	Deposit PBS scripts in each step folder for the corresponding step and make master scripts.
	'''
	
	if single == None: raise Exception('except: unclear target')
	#---single procedure script maker
	else:
		target = single
		print helpstring
		print "\tGenerating singleton script.\n"
	
		#---update the input specs file simscale parameter according to the target simulation
		simscale = 'cgmd' if any([j=='cgmd' for j in target.split('-')]) else 'aamd'
		for fn in glob.glob('./inputs/input-specs*'):
			sub = subprocess.call(['sed','-i',"s/^.*simscale.*$/simscale = \'"+simscale+"\'/",fn])
	
		#---determine location
		hostname = None
		for key in valid_hostnames.keys():
			subproc_failed = False
			try: check_host = re.match(key,subprocess.check_output(['echo $HOSTNAME'],
				shell=True).strip('\n'))
			except: check_host = False
			if check_host or re.match(key,socket.gethostname()):
				hostname = key
				break
		print '\thostname = '+str(hostname)
		#---if no hostname matches use the local steps
		if hostname == None: system_id = 'local'
		#---check for multiple architectures and choose the first one by default
		if hostname in valid_hostnames.keys() and valid_hostnames[hostname] != None: 
			arch = valid_hostnames[hostname][0]
		else: arch = None
		if hostname != None: system_id = hostname+('' if arch == None else '_'+arch)
		if hostname != None and system_id in default_proc_specs.keys():
			print '\thost/architecture settings: '
			for key in default_proc_specs[system_id]: 
				print '\t\t'+key+' = '+str(default_proc_specs[system_id][key])
			proc_settings = default_proc_specs[system_id]
			print '\tto use different settings, edit amx-bearings and rerun make'	
		else: proc_settings = None

		#---write gmxpaths.conf
		#---note that you can only change the NPROCS commands passed downstream here
		#---...this means that ./controller make must be rerun to switch processor counts
		fp = open('gmxpaths.conf','w')
		for key in standard_gromacs_commands:
			if system_id in gmx_overrides.keys() and key in gmx_overrides[system_id].keys():
				command_syntax = gmx_overrides[system_id][key]
				if proc_settings != None:
					command_syntax = re.sub('NPROCS',
						str(proc_settings['nodes']*proc_settings['ppn']),command_syntax)
				fp.write(key+' '+command_syntax+'\n')
			elif system_id in gmx_suffixes.keys(): 
				fp.write(key+' '+key+gmx_suffixes[system_id]+'\n')
			else: fp.write(key+' '+key+'\n')
		#---extra tool paths
		for key in tool_paths.keys():
			fp.write(key+' '+tool_paths[key]+'\n')
		fp.close()

		#---write the local script
		print '\twriting script: script-master-'+str(target)
		fp = open('script-master-'+str(target),'w')
		fp.write('#!/bin/bash\n')
		fp.write(script_maker(target,script_dict))
		fp.close()
		os.system('chmod u+x '+'script-master-'+str(target))
		
		#---prepare the files and directories
		prep_scripts(target,script_dict)
		
		#---for simulations that can be continued we write a script for only the final step in the sequence
		if 'continue' in script_dict[target].keys() and script_dict[target]['continue']:
			#---use chain_steps to find the final folder to deposit a continuation script if necessary
			startstep,oldsteps = chain_steps()
			fp = open(oldsteps[-1]+'/script-md-continue','w')
			fp.write('#!/bin/bash\n')
			fp.write(script_maker(target,script_dict,sim_only=True))
			fp.close()
			os.system('chmod u+x '+oldsteps[-1]+'/script-md-continue')
	
		#---write cluster-header if possible
		scratch_suffix = ''
		if proc_settings != None:
			if proc_settings['scratch']: scratch_suffix = '_scratch'
		if hostname != None and 'cluster_header_'+system_id+scratch_suffix in script_repo:
			scripttext = script_repo['cluster_header_'+system_id+scratch_suffix]
			#---if the header is a list, then it must contain a header and a footer
			if type(scripttext) == list:
				header_source = scripttext[0]
				header_source_footer = scripttext[1]
			else:
				header_source = scripttext
				header_source_footer = None
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
			print '\twriting PBS script: cluster-master-'+str(target)
			fp = open('cluster-master-'+str(target),'w')
			fp.write(header_source_mod)
			fp.write(script_maker(target,script_dict,module_commands=proc_settings['module']))
			if header_source_footer != None: fp.write(header_source_footer)
			fp.close()
			#---write a simulation continuation script to the final step
			fp = open(script_dict[target]['steps']['step_simulation']+'/cluster-md-continue','w')
			fp.write(header_source_mod)
			fp.write(script_maker(target,script_dict,module_commands=proc_settings['module'],sim_only=True))
			fp.close()		

		print '\texecute locally with ./'+'script-master-'+str(target)
		print '\tsee the documentation for details\n'

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
	argdict = {
		'clean':{'args':['protected','sure'],'module_name':None,'defaults':{'sure':False}},
		'script':{'module_name':None,'defaults':[],'args':[],
			'singles':[
				'aamd-bilayer',
				'cgmd-bilayer',
				'aamd-protein',
				'cgmd-protein',
				'protein-homology',
				'restart',
				]},
		}
	
	#---make a copy of the dictionary for pruning
	if func not in argdict.keys():
		globals()[func]()
		return
	else: argd = copy.deepcopy(argdict[func])
	
	#---special keywords are handled below
	if func == 'clean':
		for default in argd['defaults'].keys():
			kwargs[default] = argd['defaults'][default]
		for a in list(arglist): 
			kwargs[a] = True if a in argd['args'] else False
			if a in argd['args']: arglist.remove(a)
	elif func == 'script':
		if any([i for i in arglist if i in argd['singles']]) and len(arglist) > 1:
			raise Exception('except: too many arguments')
		elif any([i for i in arglist if i in argd['singles']]):
			kwargs['single'] = arglist[0]
			arglist = []

	#---all remaining keywords are handled as flags
	for a in list(arglist): 
		kwargs[a] = True if a in argd['args'] else False
		if a in argd['args']: arglist.remove(a)
	if arglist != []: raise Exception('except: unprocessed arguments: '+str(arglist))

	#---call the function possibly in another module
	if argd['module_name'] == None and func != 'gitpush': target = globals()[func]
	else: target = getattr(getattr(controller,argd['module_name']),func)
	target(**kwargs)
	return	

#---MAIN
#-------------------------------------------------------------------------------------------------------------

if __name__ == "__main__": makeface(sys.argv[1:])
	
