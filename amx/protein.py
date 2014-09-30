#!/usr/bin/python

#---AUTOMACS (AMX) "BILAYER" MODULE
#---created 2014.09.23 by ryan bradley

#---development features
import os
if os.path.isfile('/etc/pythonstart'):
	execfile('/etc/pythonstart')

#---imports
import os,sys,time,datetime,re
from tools import call,checkout,tee,copy
import amxsim

#---locate gmxpaths
if os.path.isfile('gmxpaths.conf'): gmxpaths_path = './'
else: gmxpaths_path = '../'

#---load the gromacs paths
gmxpaths = dict([[i.split()[0],' '.join(i.split()[1:])] 
	for i in open(gmxpaths_path+'gmxpaths.conf','r') 
	if i.strip()[0] != '#'])

class ProteinWater(amxsim.AMXSimulation):

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
			copy(self.sources_dir+'aamd-protein-structs/'+self.params['input_filename'],
				self.rootdir+'system-input.pdb')
			#---copy input files from standard sources i.e. the forcefield only if absent
			if needs_file_transfers:
				#---scale-specific copy commands
				if self.simscale == 'aamd':
					copy(self.sources_dir+'aamd-protein-construct/*',self.rootdir)
				else: raise Exception('except: unclear simulation resolution')
			#---running list of composition
			#---note that this method is only set to simulate a single protein
			self.nprots = 1
			self.npoz,self.nneg = 0,0
			self.nsol = 0
			#---call the master construction procedure
			print 'starting protein + water construction'
			self.construction()
			
	def construction(self):
		if not os.path.isfile(self.rootdir+'vacuum-minimized.gro'): self.vacuum()
		else: print 'skipping vacuum construction because vacuum-minimized.gro exists'

	def vacuum(self):

		#---fix histidine naming according to the convention set by the force field
		if self.settings['force_field'] == 'charmm27':
			if self.settings['histype'] == 'd':
				hisfix = "awk '{gsub(/HIS /,\"HISD\");print}' < system-input.pdb > prep-protein-start.pdb"
			if self.settings['histype'] == 'e':
				hisfix = "awk '{gsub(/HIS /,\"HISE\");print}' < system-input.pdb > prep-protein-start.pdb"
			if self.settings['histype'] == 'p':
				hisfix = "awk '{gsub(/HIS /,\"HISP\");print}' < system-input.pdb > prep-protein-start.pdb"
			call(hisfix,cwd=self.rootdir)

		print "running pdb2gmx"
		cmd = [gmxpaths['pdb2gmx'],
			'-f prep-protein-start.pdb',
			'-o vacuum-alone.gro',
			'-p vacuum-standard.top',
			'-ignh',
			'-i system-posre.itp',
			'-ff '+self.settings['force_field'],
			'-water '+self.settings['water_model']]
		call(cmd,logfile='log-pdb2gmx',cwd=self.rootdir)
		
		#---intervening step will isolate the ITP data from the TOP file to use standardized TOP
		with open(self.rootdir+'vacuum-standard.top','r') as f: topfile = f.read()	
		fp = open(self.rootdir+'protein.itp','w') 
		for line in topfile.split('\n'):
			#---skip any part of the top that follows the water topology and/or system composition
			if re.match('; Include water topology',line): break
			if re.match('; Include topology for ions',line): break
			if re.match('\[ system \]',line): break
			#---you must extract forcefield.itp from the file to prevent redundant includes
			if not re.match(".+forcefield\.itp",line) and not \
				re.match("; Include forcefield parameters",line): 
				fp.write(line+'\n')
		fp.close()
		
		#---note that this method is currently set to only simulate one protein
		self.nprots = 1
		self.write_topology_protein('vacuum.top')
		
		print "building box with "+str(self.settings['wbuffer'])+'nm of water'
		cmd = [gmxpaths['editconf'],
			'-f vacuum-alone.gro',
			'-d '+str(self.settings['wbuffer']),
			'-o vacuum.gro']
		call(cmd,logfile='log-editconf-vacuum',cwd=self.rootdir)
		
		self.minimization_method('vacuum')

		print "checking the size of the protein"
		cmd = [gmxpaths['editconf'],
			'-f vacuum-minimized.gro',
			'-o solvate-box-alone.gro',
			'-d 0']
		call(cmd,logfile='log-editconf-checksize',cwd=self.rootdir)
		boxdims = self.get_box_vectors('log-editconf-checksize')
		
		boxvecs = [i+2*self.settings['wbuffer'] for i in boxdims]
		if self.settings['cube']: 
			print "the box will be a cube"
			boxvecs = [max(boxvecs) for i in range(3)]
		print "box vectors = "+str(boxvecs)

		print "centering the protein in a new box"
		center = [i/2. for i in boxvecs]
		cmd = [gmxpaths['editconf'],
			'-f vacuum-minimized.gro',
			'-o solvate-protein.gro',
			'-center '+str(center[0])+' '+str(center[1])+' '+str(center[2]),
			'-box '+str(boxvecs[0])+' '+str(boxvecs[1])+' '+str(boxvecs[2])]
		call(cmd,logfile='log-editconf-center-protein',cwd=self.rootdir)

		print "solvating the protein in a water box"		
		copy(self.rootdir+'vacuum.top',self.rootdir+'solvate-standard.top')
		cmd = [gmxpaths['genbox'],
			'-cp solvate-protein.gro',
			'-cs spc216.gro',
			'-o solvate-dense.gro',
			'-p solvate-standard.top']
		call(cmd,logfile='log-genbox-solvate',cwd=self.rootdir)

		#---note script-construct.sh included a VMD water trimming step here
		copy(self.rootdir+'solvate-dense.gro',self.rootdir+'solvate.gro')

		cmd = ['echo -e "q\n" |',
			gmxpaths['make_ndx'],
			'-f solvate.gro',
			'-o solvate-water-check.ndx']
		call(cmd,logfile='log-make-ndx-solvate-check',cwd=self.rootdir)	
		self.nsol = int(checkout(["awk","'/ "+self.settings['sol_name']+" / {print $4}'",
			"log-make-ndx-solvate-check"],cwd=self.rootdir))/3
		self.write_topology_protein('solvate.top')

		print "minimizing with solvent"
		self.minimization_method('solvate')
		
		print "adding counterions"
		copy(self.rootdir+'solvate.top',self.rootdir+'counterions.top')
		cmd = [gmxpaths['grompp'],
			'-f input-em-steep-in.mdp',
			'-po genion.mdp',
			'-c solvate-minimized.gro',
			'-p counterions.top',
			'-o genion.tpr']
		call(cmd,logfile='log-grompp-genion',cwd=self.rootdir)
		cmd = ['echo -e "keep 0\nr '+self.settings['sol_name']+'\nkeep 1\nq\n" |',
			gmxpaths['make_ndx'],
			'-f solvate-minimized.gro',
			'-o solvate-waters.ndx']
		call(cmd,logfile='log-make-ndx-counterions-check',cwd=self.rootdir)
		cmd = [gmxpaths['genion'],
			'-s genion.tpr',
			'-o counterions.gro',
			'-conc '+str(self.settings['ion_strength']),
			'-n solvate-waters.ndx',
			'-neutral']
		call(cmd,logfile='log-genion',cwd=self.rootdir)
	
		print "updating topology with counterions"
		self.npoz,self.nneg = [int(checkout(["awk","'/Will try to add / {print $"+str(col)+\
			"}'","log-genion"],cwd=self.rootdir)) for col in [5,9]]
		self.nsol = self.nsol - self.npoz - self.nneg
		self.write_topology_protein('counterions.top')
		
		print "minimizing with solvent"
		self.minimization_method('counterions')
		
		print "completed minimization"
		copy(self.rootdir+'counterions-minimized.gro',self.rootdir+'system.gro')
		copy(self.rootdir+'counterions.top',self.rootdir+'system.top')
		self.grouping(grouptype='standard')
		

