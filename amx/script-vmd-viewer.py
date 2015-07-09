#!/usr/bin/python

"""
Quickly visualize simulations in VMD using standard representations.
"""

def uniq(a): 
	if len(a) == 1: return a[0]
	else: raise Exception('non-unique list with length '+str(len(a)))

import sys,re
import tempfile,shutil
from tools import call

#---valid coordinate/trajectory file extensions
groexts = ['gro','xtc']

#---a list of relatively self-contained functions
tcllib = {
	'orthographic':'display projection Orthographic',
	'points':['mol delrep 0 top','mol representation Points 10','mol addrep 0'],
	'protein':['mol modselect 0 top "not water"','mol modstyle 0 top NewCartoon'],
	'sideview':['mouse stoprotation','rotate x to -90','rotate y by -90'],
	'box':'pbc box',
	}

scripts = {
	'cgmd-standard':['points','orthographic','sideview','box'],
	'aamd-standard':['protein','orthographic','sideview','box'],
	}

#---default script selection
script_selection = 'cgmd-standard'

#---process arguments
args = list([i for i in sys.argv[1:] if not any([re.match('.+\.'+g+'$',i) for g in groexts])])
files = list([i for i in sys.argv[1:] if any([re.match('.+\.'+g+'$',i) for g in groexts])])
if 'cgmd' in args:
	cgmd = True
	args.remove('cgmd')
else: cgmd = False
if 'aamd' in args:
	aamd = True
	args.remove('aamd')
	script_selection = 'aamd-standard'
else: aamd = False

#---get step
if any([re.match('step=(.+)',i) for i in args]):
	stepper = int(uniq([re.findall('step=(.+)',i)[0] for i in args if re.match('step=(.+)',i)]))
else: stepper = None

#---override default script here if possible

#---assume remaining arguments are files to load
if files == []: raise Exception('no files to load')

#---script goes in temporary directory to avoid clutter
tmpdir = tempfile.mkdtemp()
print tmpdir
with open(tmpdir+'/script-vmd-look.tcl','w') as fp:
	fp.write('mol new '+files[0]+'\n')
	for fn in files[1:]:
		fp.write('mol addfile '+fn+(' step '+str(stepper) if stepper != None else '')+'\n')
	script = scripts[script_selection]
	for item in script:
		func = tcllib[item]
		if type(func) == list:
			for line in func: fp.write(line+'\n')
		elif type(func) == str: fp.write(func+'\n')
		else: raise Exception('incorrect tcllib function')

#---run the script
call('vmd -e '+tmpdir+'/script-vmd-look.tcl')
shutil.rmtree(tmpdir)
