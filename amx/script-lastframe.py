#!/usr/bin/python

#---prepare to load header
import os,sys
if os.path.abspath(os.path.curdir).split('/')[-1] != 'amx':
	#---if we execute from the root automacs directory, add amx to the path
	if os.path.isdir(os.path.abspath('./amx')): 
		prefix = 'amx/'
		sys.path.insert(0, os.path.abspath('./amx'))
	else: prefix = ''

#---load modules
execfile(prefix+'header')

#---locate gmxpaths
if os.path.isfile('gmxpaths.conf'): gmxpaths_path = './'
else: gmxpaths_path = '../'

#---load the gromacs paths
gmxpaths = dict([[i.split()[0],' '.join(i.split()[1:])] 
	for i in open(gmxpaths_path+'gmxpaths.conf','r') 
	if i.strip()[0] != '#'])

print \
"""LASTFRAME UTILITY
This program will retreive the last frame of a simulation that crashed. The user provides a file prefix"""+\
""" and the program checks for a prefix.gro file. If this file is missing, the program takes the """+\
"""prefix.xtc file and extracs the last frame. A number of log files will be written. The user must also"""+\
""" supply a root directory."""

rootdir = raw_input('root directory = ? ')
prefix = raw_input('file prefix = ? ')
lastframe(rootdir=rootdir,prefix=prefix,gmxpaths=gmxpaths)

