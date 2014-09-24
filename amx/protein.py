#!/usr/bin/python

#---AUTOMACS (AMX) "BILAYER" MODULE
#---created 2014.09.23 by ryan bradley

#---development features
import os
if os.path.isfile('/etc/pythonstart'):
	execfile('/etc/pythonstart')

#---imports
import os
import sys
import time
import datetime
from tools import call,checkout,tee,copy

#---locate gmxpaths
if os.path.isfile('gmxpaths.conf'): gmxpaths_path = './'
else: gmxpaths_path = '../'

#---load the gromacs paths
gmxpaths = dict([[i.split()[0],' '.join(i.split()[1:])] 
	for i in open(gmxpaths_path+'gmxpaths.conf','r') 
	if i.strip()[0] != '#'])

class ProteinWater:

	def __init__(self,rootdir=None,**kwargs):
		#---set paths for the previous and current steps
		self.rootdir = os.path.abspath(os.path.expanduser(rootdir))+'/'

		#---manually specify the sources directory and the inputs file
		if 'sources_dir' in kwargs.keys(): 
			self.sources_dir = os.path.abspath(os.path.expanduser(kwargs['sources_dir']))+'/'
		else: self.sources_dir = os.path.abspath(os.path.expanduser('./sources/'))+'/'
		if 'inputs_file' in kwargs.keys(): 
			self.inputs_file = os.path.abspath(os.path.expanduser(kwargs['inputs_file']))+'/'
		else: self.inputs_file = os.path.abspath(os.path.expanduser('./inputs/input-specs-protein.dat'))

		#---skip everything if the procedure is final configuration is already available
		#---note that this takes the place of an exception and hence expects the user to remove old files
		if os.path.isfile(self.rootdir+'system.gro'): 
			sys.stdout = tee(open(self.rootdir+'log-script-master','a',1))
			sys.stderr = tee(open(self.rootdir+'log-script-master','a',1),error=True)
			print 'skipping this step because system.gro exists'
		else:
			#---make root directory
			if not os.path.isdir(self.rootdir): 
				os.mkdir(rootdir)
				needs_file_transfers = True
			else: needs_file_transfers = False
			#---start the logger
			sys.stdout = tee(open(self.rootdir+'log-script-master','a',1))
			sys.stderr = tee(open(self.rootdir+'log-script-master','a',1),error=True)
			print 'START PROTEIN + WATER'
			#---protein construction settings are pulled from the inputs_file
			self.params = dict()
			execfile(self.inputs_file,self.params)
			self.simscale = self.params['simscale']
			self.settings = self.params['protein_construction_settings'][self.simscale]
			#---copy input files if this is a fresh start
			copy(self.sources_dir+'aamd-protein/'+self.params['input_filename'],
				self.rootdir+'system-input.pdb')
			#---copy input files from standard sources i.e. the forcefield only if absent
			if needs_file_transfers:
				#---scale-specific copy commands
				if self.simscale == 'aamd':
					copy(self.sources_dir+'aamd-protein-construct/*',self.rootdir)
				else: raise Exception('except: unclear simulation resolution')
			#---running lists of molecules and compositions
			self.lnames,self.comps = [],[]
			#---call the master construction procedure
			print 'starting protein + water construction'
			self.construction()
			
	def minimization_method(self,name,posre=False):
		'''
		Generic minimization method with a drop-in name and standard inputs. This function takes a file
		suffix and performs a minimization on the corresponding ``gro`` file.
		'''
		print name+" minimization, steepest descent"
		cmd = [gmxpaths['grompp'],
			'-f input-em-steep-'+('posre-' if posre else '')+'in.mdp',
			'-c '+name+'.gro',
			'-p '+name+'.top',
			'-o em-'+name+'-steep',
			'-po em-'+name+'-steep',
			'-maxwarn 10']
		call(cmd,logfile='log-grompp-em-'+name+'-steep',cwd=self.rootdir)
		cmd = [gmxpaths['mdrun'],'-v','-deffnm em-'+name+'-steep']
		call(cmd,logfile='log-mdrun-em-'+name+'-steep',cwd=self.rootdir)

		print name+" minimization, conjugate gradient"
		cmd = [gmxpaths['grompp'],
			'-f input-em-cg-'+('posre-' if posre else '')+'in.mdp',
			'-c em-'+name+'-steep.gro',
			'-p '+name+'.top',
			'-o em-'+name+'-cg',
			'-po em-'+name+'-cg',
			'-maxwarn 10']
		print ' '.join(cmd)
		call(cmd,logfile='log-grompp-em-'+name+'-cg',cwd=self.rootdir)
		cmd = [gmxpaths['mdrun'],'-v','-deffnm em-'+name+'-cg']
		call(cmd,logfile='log-mdrun-em-'+name+'-cg',cwd=self.rootdir)
	
		mdrun_logs = [[line for line in open(self.rootdir+fn,'r')] for fn in 
			['log-mdrun-em-'+name+'-steep','log-mdrun-em-'+name+'-cg']]
		forces = [i[0] if len(i) >= 1 else [] for i in [[float(line.strip().split()[3]) 
			for line in d if line[:13] == 'Maximum force'] for d in mdrun_logs]]
		if len(forces) == 0: raise Exception('except: both minimization methods failed')
		if forces[0] == []: bestmin = 'em-'+name+'-cg.gro'
		elif forces[1] == []: bestmin = 'em-'+name+'-steep.gro'
		else: bestmin = 'em-'+name+'-cg.gro' if forces[0] > forces[1] else 'em-'+name+'-steep.gro'
		call('cp '+bestmin+' '+name+'-minimized.gro',cwd=self.rootdir)
	
	def construction(self):
		if not os.path.isfile(self.rootdir+'vacuum-minimized.gro'): self.vacuum()
		else: print 'skipping vacuum construction because vacuum-minimized.gro exists'

	def vacuum(self):

		#---fix histidine naming according to the convention set by the force field
		if self.settings['force_field'] == 'charmm27':
			if self.settings['histype'] == 'd':
				hisfix = "awk '{gsub(/HIS /,\"HID \");print}' < system-input.pdb > prep-protein-start.pdb"
			if self.settings['histype'] == 'e':
				hisfix = "awk '{gsub(/HIS /,\"HSE \");print}' < system-input.pdb > prep-protein-start.pdb"
			if self.settings['histype'] == 'p':
				hisfix = "awk '{gsub(/HIS /,\"HSP \");print}' < system-input.pdb > prep-protein-start.pdb"
			call(hisfix,cwd=self.rootdir)

		print "running pdb2gmx"
		cmd = [gmxpaths['pdb2gmx'],
			'-f prep-protein-start.pdb',
			'-o vacuum-alone.gro',
			'-p vacuum.top',
			'-ignh',
			'-i system-posre.itp',
			'-ff '+self.settings['force_field'],
			'-water '+self.settings['water_model'],
			]
		print ' '.join(cmd)
		call(cmd,logfile='log-pdb2gmx',cwd=self.rootdir)

		print "building box with "+str(self.settings['wbuffer'])+'nm of water'
		cmd = [gmxpaths['editconf'],
			'-f vacuum-alone.gro',
			'-d '+str(self.settings['wbuffer']),
			'-o vacuum.gro',
			]
		print ' '.join(cmd)
		call(cmd,logfile='log-editconf-vacuum',cwd=self.rootdir)
		
		self.minimization_method('vacuum')
		
		#--->>>> UNDER DEVELOPMENT AFTER HERE
		
		print "running "
		cmd = [gmxpaths[''],
			]
		print ' '.join(cmd)
		call(cmd,logfile='log-',cwd=self.rootdir)
		
		
		
		
		
		
		
