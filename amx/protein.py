#!/usr/bin/python

#---AUTOMACS (AMX) "BILAYER" MODULE
#---created 2014.09.23 by ryan bradley

#---development features
import os,subprocess
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
			#---the presence of the root directory signals that the files have already been copied
			#---...note that we prefer to allow the user to modify input files and restart so it is
			#---...important to retain this feature which is why we don't pre-make these directories
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
			if os.path.isfile(self.sources_dir+'aamd-protein-structs/'+self.params['input_filename']):
				copy(self.sources_dir+'aamd-protein-structs/'+self.params['input_filename'],
					self.rootdir+'system-input.pdb')
			elif os.path.isfile(os.path.abspath(self.params['input_filename'])):
				copy(self.params['input_filename'],self.rootdir+'system-input.pdb')
			else: raise Exception('except: cannot locate input_filename '+self.params['input_filename'])
			#---copy input files from standard sources i.e. the forcefield only if absent
			if needs_file_transfers:
				#---scale-specific copy commands
				if self.simscale == 'aamd':
					copy(self.sources_dir+'aamd-protein-construct/*',self.rootdir)
				elif self.simscale == 'cgmd':
					copy(self.sources_dir+'cgmd-protein-construct/*',self.rootdir)
					copy(self.sources_dir+'martini.ff',self.rootdir+'martini.ff')
				else: raise Exception('except: unclear simulation resolution')
			#---running list of composition
			#---note that this method is only set to simulate a single protein
			self.nprots = 1
			self.npoz,self.nneg = 0,0
			self.nsol = 0
			self.protname = 'Protein'
			#---call the master construction procedure
			print 'starting protein + water construction'
			self.construction()
			
	def construction(self):
		if self.simscale == 'aamd':
			if not os.path.isfile(self.rootdir+'vacuum-minimized.gro'): self.vacuum()
			else: print 'skipping vacuum construction because vacuum-minimized.gro exists'
		elif self.simscale == 'cgmd':	
			if not os.path.isfile(self.rootdir+'vacuum-minimized.gro'): self.vacuum_cgmd()
			else: print 'skipping vacuum construction because vacuum-minimized.gro exists'
		else: raise Exception('except: unclear simulation scale')
		if not os.path.isfile(self.rootdir+'solvate-minimized.gro'): self.solvate()
		else: print 'skipping solvate construction because solvate-minimized.gro exists'
		if not os.path.isfile(self.rootdir+'counterions-minimized.gro'): self.counterions()
		else: print 'skipping counterions construction because counterions-minimized.gro exists'
		if not os.path.isfile(self.rootdir+'system.gro'): self.groups()
		else: print 'skipping the grouping step because system.gro exists'
	
	def vacuum_cgmd(self):

		exstring_dssp = 'except: cannot find dssp at '+gmxpaths['dssp']+\
			'\nconsider using the following syntax to download for 64-bit linux:'+\
			'\n\twget ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.0.4-linux-amd64'+\
			'\n\tor navigate to ftp://ftp.cmbi.ru.nl/pub/software/dssp/'
			
		exstring_martinize = 'except: cannot find martinize at '+gmxpaths['martinize']+\
			'\nconsider using the following syntax to download:'+\
			'\n\twget http://md.chem.rug.nl/cgmartini/images/tools/martinize/martinize-2.4/martinize.py'+\
			'\n\tor navigate to http://md.chem.rug.nl/cgmartini/index.php/tools2/proteins-and-bilayers'
	
		#---first test to see if executables are available
		if not os.path.isfile(os.path.expanduser(gmxpaths['dssp'])): raise Exception(exstring_dssp)
		if not os.path.isfile(os.path.expanduser(gmxpaths['martinize'])): raise Exception(exstring_martinize)	
	
		cmd = [gmxpaths['martinize'],
			'-f system-input.pdb',
			'-o system-original.top',
			'-x protein-cg.pdb',
			'-ff martini22','-ed',
			'-dssp '+gmxpaths['dssp']]
		call(cmd,logfile='log-martinize',cwd=self.rootdir)
		
		with open(self.rootdir+'system-original.top') as fp: lines = fp.readlines()
		self.itp_protein = [l.split()[0] for l in lines if l[:7] == 'Protein']

		#---note that this section leaves out lipids
		self.itp_lipid = []
		
		#---note that this method is currently set to only simulate one protein
		self.nprots = [1]
		self.write_topology_protein('vacuum.top')
		
		cmd = [gmxpaths['editconf'],
			'-f protein-cg.pdb',
			'-o vacuum-alone.gro']
		call(cmd,logfile='log-editconf-convert',cwd=self.rootdir)
	
		print "building box with "+str(self.settings['wbuffer'])+'nm of water'
		cmd = [gmxpaths['editconf'],
			'-f vacuum-alone.gro',
			'-d '+str(self.settings['wbuffer']),
			'-o vacuum.gro','-c']
		call(cmd,logfile='log-editconf-vacuum',cwd=self.rootdir)
		
		self.minimization_method('vacuum')

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
			
		print "stripping non-protein molecules"
		cmd = [gmxpaths['make_ndx'],
			'-f prep-protein-start.pdb',
			'-o prep-index-protein.ndx']
		call(cmd,logfile='log-make-ndx-prep-protein',cwd=self.rootdir,inpipe="q\n")
		protgrp = int(checkout(["awk","'/[ ,\t]+Protein[ ,\t]+:/ {print $1}'",
			"log-make-ndx-prep-protein"],cwd=self.rootdir).strip())
		cmd = [gmxpaths['make_ndx'],
			'-f prep-protein-start.pdb',
			'-o prep-index-protein-only.ndx']
		call(cmd,logfile='log-make-ndx-prep-protein-only',cwd=self.rootdir,
			inpipe="keep "+str(protgrp)+"\nq\n")
		cmd = [gmxpaths['editconf'],
			'-f prep-protein-start.pdb',
			'-o prep-protein-start-stripped.pdb',
			'-n prep-index-protein-only.ndx']
		call(cmd,logfile='log-editconf-prep-protein-strip',cwd=self.rootdir)

		print "running pdb2gmx"
		cmd = [gmxpaths['pdb2gmx'],
			'-f prep-protein-start-stripped.pdb',
			'-o vacuum-alone.gro',
			'-p vacuum-standard.top',
			'-ignh',
			'-i system-posre.itp',
			'-ff '+self.settings['force_field'],
			'-water '+self.settings['water_model']]
		call(cmd,logfile='log-pdb2gmx',cwd=self.rootdir)
		
		#---intervening step will isolate the ITP data from the TOP file to use standardized TOP
		with open(self.rootdir+'vacuum-standard.top','r') as f: topfile = f.read()
		#---extract protein chain names here if necessary
		chains = []
		startline = [ii for ii,i in enumerate(topfile.split('\n')) 
			if re.match('^(\s+)?\[(\s+)?system(\s+)?\]',i)][0]
		for line in topfile.split('\n')[startline:]:
			if re.match('^Protein',line):
				chains.append(line.split(' ')[0])
		if len(chains) > 1:
			#---assume one domain per chain
			self.nprots = [1 for i in chains]
			self.protname = chains
		else:	
			self.protname = chains[0]
			self.nprots = 1
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
		self.write_topology_protein('vacuum.top')
		
		print "building box with "+str(self.settings['wbuffer'])+'nm of water'
		cmd = [gmxpaths['editconf'],
			'-f vacuum-alone.gro',
			'-d '+str(self.settings['wbuffer']),
			('-princ' if 'align_x' in self.settings.keys() 
			and self.settings['align_x'] == True else ''),
			'-o vacuum.gro']
		call(cmd,logfile='log-editconf-vacuum',cwd=self.rootdir,
			inpipe=('0\n' if 'align_x' in self.settings.keys() 
			and self.settings['align_x'] == True else None))		
		self.minimization_method('vacuum')

	def solvate(self):

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
		
		if self.simscale == 'aamd': 

			#---note script-construct.sh included a VMD water trimming step here
			print "solvating the protein in a water box"		
			copy(self.rootdir+'vacuum.top',self.rootdir+'solvate-standard.top')
			cmd = [gmxpaths['genbox'],
				'-cp solvate-protein.gro',
				'-cs '+self.settings['solvent_structure'],
				'-o solvate-dense.gro',
				'-p solvate-standard.top']
			call(cmd,logfile='log-genbox-solvate',cwd=self.rootdir)
			copy(self.rootdir+'solvate-dense.gro',self.rootdir+'solvate.gro')

		elif self.simscale == 'cgmd':

			#---manually create a large water box
			print 'creating a custom water box'
			nmult = [int(i/3.644+1) for i in boxvecs]
			cmd = [gmxpaths['genconf'],
				'-f '+self.settings['solvent_structure'],
				'-o solvate-water-big.gro',
				'-nbox '+' '.join([str(i) for i in nmult])]
			call(cmd,logfile='log-genconf-replicate',cwd=self.rootdir)
			#---concatenate the protein and water box
			with open(self.rootdir+'solvate-protein.gro','r') as fp: gro1 = fp.readlines()
			with open(self.rootdir+'solvate-water-big.gro','r') as fp: gro2 = fp.readlines()
			natoms = int(gro1[1].strip())+int(gro2[1].strip())
			with open(self.rootdir+'solvate-merge.gro','w') as fp:
				fp.write('protein+water merged\n')
				fp.write(str(natoms)+'\n')
				for line in gro1[2:-1]: fp.write(line)
				for line in gro2[2:]: fp.write(line)		
			vmdtrim = [
				'mol new solvate-merge.gro',
				'set sel [atomselect top \"(all not (name '+self.settings['sol_name']+\
				' and within '+str(self.settings['protein_water_gap'])+\
				' of not name '+self.settings['sol_name']+')) and '+\
				'(x>=0 and x<='+str(10*boxvecs[0])+' and y>=0 and y<= '+str(10*boxvecs[1])+\
				' and z>=0 and z<= '+str(10*boxvecs[2])+')"]',
				'$sel writepdb '+self.rootdir+'solvate-vmd.pdb',
				'exit',]			
			with open(self.rootdir+'script-vmd-trim.tcl','w') as fp:
				for line in vmdtrim: fp.write(line+'\n')
			vmdlog = open(self.rootdir+'log-script-vmd-trim','w')
			p = subprocess.Popen('vmd -dispdev text -e script-vmd-trim.tcl',
				stdout=vmdlog,stderr=vmdlog,cwd=self.rootdir,shell=True)
			p.communicate()
			cmd = [gmxpaths['editconf'],
				'-f solvate-vmd.pdb',
				'-o solvate.gro','-resnr 1','-d 0']
			call(cmd,logfile='log-editconf-convert',cwd=self.rootdir)

		cmd = [gmxpaths['make_ndx'],
			'-f solvate.gro',
			'-o solvate-water-check.ndx']
		call(cmd,logfile='log-make-ndx-solvate-check',cwd=self.rootdir,inpipe="q\n")	
		self.nsol = int(checkout(["awk","'/ "+self.settings['sol_name']+" / {print $4}'",
			"log-make-ndx-solvate-check"],cwd=self.rootdir))/(3 if self.simscale == 'aamd' else 1)
		self.write_topology_protein('solvate.top')

		print "minimizing with solvent"
		self.minimization_method('solvate')

	def counterions(self):

		print "adding counterions"
		copy(self.rootdir+'solvate.top',self.rootdir+'counterions.top')
		cmd = [gmxpaths['grompp'],
			'-f input-em-steep-in.mdp',
			'-po genion.mdp',
			'-c solvate-minimized.gro',
			'-p counterions.top',
			'-o genion.tpr']
		call(cmd,logfile='log-grompp-genion',cwd=self.rootdir)
		cmd = [gmxpaths['make_ndx'],
			'-f solvate-minimized.gro',
			'-o solvate-waters.ndx']
		call(cmd,logfile='log-make-ndx-counterions-check',cwd=self.rootdir,
			inpipe='keep 0\nr '+self.settings['sol_name']+'\nkeep 1\nq\n')
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

	def groups(self):

		print "completed minimization"
		copy(self.rootdir+'counterions-minimized.gro',self.rootdir+'system.gro')
		copy(self.rootdir+'counterions.top',self.rootdir+'system.top')
		if self.simscale == 'aamd': grouptype = 'standard'
		if self.simscale == 'cgmd': grouptype = 'cgmd_water'
		self.grouping(grouptype=grouptype)
		

