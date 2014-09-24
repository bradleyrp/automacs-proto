#!/usr/bin/python

#---AUTOMACS (AMX) "BILAYER" MODULE
#---created 2014.08.14 by ryan bradley

"""
This module takes one or two monolayers consisting of a grid of lipids oriented in the z-dimension. It
combines them into a bilayer, minimizes in vacuum, simulates with restraints to pack the lipids into a natural
area, adds water and counterions, and performs a final minimization. The resulting structure will be ready for
equilibration and production simulation.

.. figure:: ../sources/docs/automacs-cgmd-bilayer-bilayer.jpeg
   :align: center
   :alt: monolayer snapshot
   :width: 75 %
   
   Example bilayer configuration before adding solvent.

For simplicity, you build a bilayer by creating a ``Bilayer`` class instance. You must provide the location of
the monolayer configurations necessary to build the bilayer. This can be done by passing the ``previous_dir``
keyword argument according to the following example from ``script-amx-cgmd-bilayer``. ::

	Bilayer(rootdir='s2-build-bilayer',previous_dir='s1-build-lipidgrid')
	
Creating the class object takes care of everything else. When it is complete, you will find a ``system.gro`` 
and the associated topology and index files in the ``s2-build-bilayer`` folder. The assembly process includes
the following steps.

1. Assemble the monolayers into a bilayer and minimize.
2. Pack the lipids into a more natural configuration by simulating in vacuum with position restraints
   on the head and tails to keep the lipids aligned.
3. Add a slab of water and minimize.
4. Add counterions and minimize.

The overall procedure is executed via ``Bilayer.construction()`` and can be restarted in case of errors at the 
last successful step.
"""

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

#---CLASS
#-------------------------------------------------------------------------------------------------------------

class Bilayer:
	'''
	Class which generates a solvated bilayer with counterions from monolayer configurations.
	
	Note that creating an instance of this class will automatically run the construction procedure. You must 
	provide locations to the monolayer configurations when the object is created. Otherwise, all parameters
	are specified in ``inputs/input-specs-bilayer.dat``. Note that this file handles inputs for both the 
	atomistic and coarse-grained bilayer construction, both of which are handled by this class. The 
	construction process procedes automatically, but will skip ahead if the process is interrupted and 
	then restarted.

	Parameters
	----------
	rootdir : str
		relative path to store simulation data
	previous_dir : str
		path to the directory containing monolayer configurations in either ``place-grid-start0.gro`` and 
		``place-grid-start0.gro`` for a bilayer with an asymmetric composition, or 
		``place-grid-start.gro`` for a symmetric bilayer.
	'''
	def __init__(self,rootdir=None,previous_dir=None,**kwargs):
		#---set paths for the previous and current steps
		self.rootdir = os.path.abspath(os.path.expanduser(rootdir))+'/'
		self.prevdir = os.path.abspath(os.path.expanduser(previous_dir))+'/'

		#---manually specify the sources directory and the inputs file
		if 'sources_dir' in kwargs.keys(): 
			self.sources_dir = os.path.abspath(os.path.expanduser(kwargs['sources_dir']))+'/'
		else: self.sources_dir = os.path.abspath(os.path.expanduser('./sources/'))+'/'
		if 'inputs_file' in kwargs.keys(): 
			self.inputs_file = os.path.abspath(os.path.expanduser(kwargs['inputs_file']))+'/'
		else: self.inputs_file = os.path.abspath(os.path.expanduser('./inputs/input-specs-bilayer.dat'))

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
			print 'START BILAYER'
			#---determine if asymmetric from presence of place-grid-start0.gro
			if os.path.isfile(self.prevdir+'place-grid-start0.gro'): self.asymmetric = True
			else: self.asymmetric = False
			#---bilayer construction settings are pulled from the inputs_file
			self.params = dict()
			execfile(self.inputs_file,self.params)
			self.simscale = self.params['simscale']
			self.settings = self.params['bilayer_construction_settings'][self.simscale]
			self.ion_residue_names = [self.settings['sol_name']]+\
				self.params['bilayer_construction_settings'][self.simscale]['ion_residue_names']
			#---copy input files if this is a fresh start
			if self.asymmetric:
				copy(self.prevdir+'place-grid-start0.gro',self.rootdir+'place-grid-start0.gro')
				copy(self.prevdir+'place-grid-start1.gro',self.rootdir+'place-grid-start1.gro')
			else: 
				copy(self.prevdir+'place-grid-start.gro',self.rootdir+'place-grid-start0.gro')
				copy(self.prevdir+'place-grid-start.gro',self.rootdir+'place-grid-start1.gro')
			copy(self.prevdir+'composition.dat',self.rootdir+'composition.dat')
			#---copy input files from standard sources i.e. the forcefield only if absent
			if needs_file_transfers:
				#---scale-specific copy commands
				if self.simscale == 'cgmd':			
					copy(self.sources_dir+'cgmd-bilayer-lipids-tops',self.rootdir+'lipids-tops')
					copy(self.sources_dir+'martini.ff',self.rootdir+'martini.ff')
					copy(self.sources_dir+'cgmd-bilayer-construct/*',self.rootdir)
					copy(self.sources_dir+'cgmd-bilayer/solvate-water.gro',self.rootdir+'solvate-water.gro')
				elif self.simscale == 'aamd':
					copy(self.sources_dir+'aamd-bilayer-lipids-tops',self.rootdir+'lipids-tops')
					copy(self.sources_dir+'charmm36.ff',self.rootdir+'charmm36.ff')
					copy(self.sources_dir+'aamd-bilayer-construct/*',self.rootdir)
				else: raise Exception('except: unclear simulation resolution')
			#---running lists of molecules and compositions
			self.lnames,self.comps = [],[]
			#---call the master construction procedure
			print 'starting bilayer construction'
			self.construction()

	def get_box_vectors(self,logfile,new=True):
		'''Read box vectors from editconf output stored in a log file.'''
		if new:
			boxdims = [float(checkout(["awk","'/new box vectors/ {print $"+\
				str(i+5)+"}'",logfile],cwd=self.rootdir).strip()) for i in range(3)]
		else:
			boxdims = [float(checkout(["awk","'/    box vectors/ {print $"+\
				str(i+4)+"}'",logfile],cwd=self.rootdir).strip()) for i in range(3)]
		return boxdims
		
	def write_topology(self,topname,water_only=False):
		'''Write the topology file from residue lists stored in ``self.lnames`` and ``self.complist``.'''
		
		print 'writing topology'
		#---atomistic topology disabled here
		if self.simscale == 'aamd':
			fp = open(self.rootdir+topname,'w')
			fp.write('#include "charmm36.ff/forcefield.itp"\n')
			for lname in self.lnames:
				if lname not in self.ion_residue_names: 
					fp.write('#include "./lipids-tops/lipid.'+lname+'.itp"\n')
					fp.write('#include "./lipids-tops/restr.'+lname+'.itp"\n') 
			fp.write('#include "charmm36.ff/tips3p.itp"\n')
			fp.write('#include "charmm36.ff/ions.itp"\n\n')
			fp.write('[ system ]\nBILAYER\n\n')
			fp.write('[ molecules ]\n')
			for l in range(len(self.lnames)): 
				if not water_only or self.lnames[l] == self.settings['sol_name']:
					fp.write(self.lnames[l]+' '+self.comps[l][1]+'\n')
			fp.close()
		elif self.simscale == 'cgmd':
			#---define ion names for exclusion from topology update function
			self.ion_residue_names = [self.settings['sol_name'],'MG','NA','CL','Cal']
			fp = open(self.rootdir+topname,'w')
			fp.write(self.params['topheader_martini'])			
			for l in range(len(self.lnames)): 
				if not water_only or self.lnames[l] == self.settings['sol_name']: 
					fp.write(self.lnames[l]+' '+self.comps[l][1]+'\n')
			fp.close()
		else: raise Exception('except: unclear procedure to write topology')
	
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
		'''
		Execute construction steps in sequence and skip ahead if they are already completed. Infer whether 
		construction steps are complete by checking for their final output files. The procedure executes in
		the following order.
		
		1. :class:`vacuum <Bilayer.vacuum>`
		2. :class:`packing <Bilayer.packing>`
		3. :class:`solvate <Bilayer.solvate>`
		4. :class:`counterionize <Bilayer.counterionize>`
		
		For development or tweaking the input parameters, the system checks for the following files at these
		corresponding stages.
		
		1. vacuum-minimized.gro
		2. vacuum-packed.gro
		3. solvate-minimized.gro
		4. counterions-minimized.gro
		5. system.gro
		
		If the system.gro file is absent, the script will regenerate the topology and group files. You may 
		delete any of these key files to repeat a process when you re-run the script.
		'''
		if not os.path.isfile(self.rootdir+'vacuum-minimized.gro'): self.vacuum()
		else: 
			print 'skipping vacuum construction because vacuum-minimized.gro exists'
			#---diagnostic
			self.comps = [line.strip().split() for line in open(self.rootdir+'composition.dat')]
			self.lnames = [i[0] for i in self.comps]
			monos = [[line for line in open(self.rootdir+fn,'r')] 
				for fn in ['place-grid-monolayer-center0.gro','place-grid-monolayer-mirror1.gro']]
		if not os.path.isfile(self.rootdir+'vacuum-packed.gro'): self.packing()
		else: print 'skipping vacuum packing because vacuum-packed.gro exists'
		if not os.path.isfile(self.rootdir+'solvate-minimized.gro'): self.solvate()
		else:
			print 'skipping solvate because solvate-minimized.gro exists'
			#---diagnostic
			if self.simscale == 'aamd':
				nwaters = (int(checkout(['wc','-l','solvate-empty-uncentered.gro'],
					cwd=self.rootdir).split()[0])-3)/3
			elif self.simscale == 'cgmd':
				nwaters = (int(checkout(['wc','-l','solvate-empty-uncentered.gro'],
					cwd=self.rootdir).split()[0])-3)				
			self.lnames.append(self.settings['sol_name'])
			self.comps.append([self.settings['sol_name'],str(nwaters)])
		if not os.path.isfile(self.rootdir+'counterions-minimized.gro'): self.counterionize()
		else: 
			print 'skipping counterionize because counterions-minimized.gro exists'
			#---diagnostic
			pcount,ncount = [int(checkout(['awk',"'/Will try to add / {print $"+str(col)+\
				"}'","log-genion"],cwd=self.rootdir)) for col in [5,9]]
			self.lnames.append(self.settings['positive_ion_name'])
			self.comps.append([self.settings['positive_ion_name'],str(pcount)])
			self.lnames.append(self.settings['negative_ion_name'])
			self.comps.append([self.settings['negative_ion_name'],str(ncount)])
			self.comps[self.lnames.index(self.settings['sol_name'])][1] = str(nwaters - pcount - ncount)
		if not os.path.isfile(self.rootdir+'system.gro'): self.grouping()
	
	def vacuum(self):	
		'''Assemble monolayers into a bilayer and minimize.'''
		startconfs = ['place-grid-start0.gro','place-grid-start1.gro']
		for si in range(len(startconfs)):
			startconf = startconfs[si]

			print "checking the unsized monolayer dimensions"
			call('cp '+startconf+' place-grid-monolayer-unsized'+str(si)+'.gro',cwd=self.rootdir)
			cmd = [gmxpaths['editconf'],
				'-f place-grid-monolayer-unsized'+str(si)+'.gro',
				'-o place-grid-monolayer-check-dims'+str(si)+'.gro',
				'-d 0']
			call(cmd,logfile='log-editconf-place-check-dims-monolayer'+str(si),cwd=self.rootdir)
			boxdims = self.get_box_vectors('log-editconf-place-check-dims-monolayer'+str(si))

			print "resizing the box with a buffer"
			cmd = [gmxpaths['editconf'],
				'-f place-grid-monolayer-unsized'+str(si)+'.gro',
				'-o place-grid-monolayer'+str(si)+'.gro',
				'-d '+str(self.settings['lbuffer'])]
			call(cmd,logfile='log-editconf-place-grid-monolayer'+str(si),cwd=self.rootdir)

			shifter = boxdims[2]+self.settings['tbuffer']
			cmd = [gmxpaths['editconf'],
				'-f place-grid-monolayer'+str(si)+'.gro',
				'-o place-grid-monolayer-center'+str(si)+'.gro',
				'-rotate 0 0 90' if si == 0 else '-rotate 0 0 0',
				'-center 0 0 0']
			call(cmd,logfile='log-editconf-place-center'+str(si),cwd=self.rootdir)

			print "rotating and flipping"
			cmd = [gmxpaths['editconf'],
				'-f place-grid-monolayer-center'+str(si)+'.gro',
				'-o place-grid-monolayer-mirror'+str(si)+'.gro',
				'-scale 1 1 -1',
				'-translate 0 0 -'+str(shifter)]
			call(cmd,logfile='log-editconf-place-mirror'+str(si),cwd=self.rootdir)

		#---diagnostic
		self.comps = [line.strip().split() for line in open(self.rootdir+'composition.dat')]
		self.lnames = [i[0] for i in self.comps]
		monos = [[line for line in open(self.rootdir+fn,'r')] 
			for fn in ['place-grid-monolayer-center0.gro','place-grid-monolayer-mirror1.gro']]

		fp = open(self.rootdir+'place-grid-unsized.gro','w')
		natoms = str(sum([len(mono)-3 for mono in monos]))
		fp.write('bilayer\n'+natoms+'\n')
		for lname in self.lnames:
			for mono in monos:
				for line in mono:
					if line[5:10].strip() == lname: 
						fp.write(line)
		fp.write(monos[0][-1])
		fp.close()

		print "checking the unsized bilayer dimensions"
		cmd = [gmxpaths['editconf'],
			'-f place-grid-unsized.gro',
			'-o place-grid-unsized-check-dims.gro',
			'-d '+str(self.settings['lbuffer'])]
		call(cmd,logfile='log-editconf-place-check-dims',cwd=self.rootdir)

		boxdims = self.get_box_vectors('log-editconf-place-check-dims',new=True)
		if self.settings['square']:
			boxpad = [max(boxdims[:2])+self.settings['lbuffer']/2. for i in range(2)]
		else:
			boxpad = [boxdims[i]+self.settings['lbuffer']/2. for i in range(2)]

		print "resizing the box" 
		cmd = [gmxpaths['editconf'],
			'-f place-grid-unsized.gro',
			'-o vacuum.gro',
			'-box '+str(boxpad[0])+' '+str(boxpad[1])+' '+str(self.settings['zbuffer']/2.+boxdims[2])]
		call(cmd,logfile='log-editconf-place-grid',cwd=self.rootdir)
		
		self.write_topology('vacuum.top')
		self.minimization_method('vacuum',posre=True)
		
	def packing(self):
		'''
		Pack the lipids into a natural configuration by simulating in vacuum with restraints.
		
		This function runs vacuum simulation in stages according to available input files specified in the
		``inputs/input-specs-bilayer.dat`` file. These inputs are sequenced to ensure that the bilayer is 
		contiguous and the lipids are aligned. Position restraints written into the lipid topologies prevent
		the lipids from escaping.
		'''
		#---find the list of mdp files used for the stagewise simulation	
		packing_mdps = self.settings['vacuum_packing_sequence']
			
		#---loop over stages in order to gently but firmly assemble the bilayer
		for mi in range(1,len(packing_mdps)+1):
			if os.path.isfile(self.rootdir+'md-vacuum-p'+str(mi)+'.gro'):
				print 'skipping packing step '+str(mi)+' because md-vacuum-p'+str(mi)+'.gro exists'
			else:
				print "packing lipids in vacuum, round "+str(mi)
				cmd = [gmxpaths['grompp'],
					'-f '+packing_mdps[mi-1],
					'-c '+('vacuum-minimized.gro' if mi == 1 else 'md-vacuum-p'+str(mi-1)+'.gro'),
					'-p vacuum.top',
					'-o md-vacuum-p'+str(mi),
					'-po md-vacuum-p'+str(mi),
					'-maxwarn 10']
				call(cmd,logfile='log-grompp-md-vacuum-p'+str(mi),cwd=self.rootdir)
				cmd = [gmxpaths['mdrun'],'-v','-deffnm md-vacuum-p'+str(mi)]
				call(cmd,logfile='log-mdrun-md-vacuum-p'+str(mi),cwd=self.rootdir)
		call('cp md-vacuum-p'+str(mi)+'.gro vacuum-packed.gro',cwd=self.rootdir)

	def solvate(self):
		'''
		Add a slab of water next to the bilayer. ``inputs/input-specs-bilayer.dat`` sets some of the geometric
		parameters, namely the solvent_thickness, which determines the total amount of water in the system 
		before any relaxation steps.
		'''
	
		#---set the file name of the water box
		#---note that spc216.gro used for atomistic water is found the GROMACS share directory
		water_conf = self.settings['water_conf']
	
		print "checking the bilayer dimensions"
		cmd = [gmxpaths['editconf'],
			'-f vacuum-packed.gro',
			'-o md-vacuum-p2-check-dims.gro',
			'-d 0']
		call(cmd,logfile='log-editconf-vacuum-packed-check-dims',cwd=self.rootdir)
		boxdims = self.get_box_vectors('log-editconf-vacuum-packed-check-dims',new=False)

		print "making a water box with the same footprint as the bilayer"
		cmd = [gmxpaths['genbox'],
			'-cs '+self.settings['water_conf'],
			'-o solvate-empty-uncentered.gro',
			'-box '+str(boxdims[0])+' '+str(boxdims[1])+' '+str(self.settings['solvent_thickness'])]
		call(cmd,logfile='log-genbox-solvate-empty',cwd=self.rootdir)

		print "counting waters"
		if self.simscale == 'aamd':
			nwaters = (int(checkout(['wc','-l','solvate-empty-uncentered.gro'],
				cwd=self.rootdir).split()[0])-3)/3
		elif self.simscale == 'cgmd':
			nwaters = (int(checkout(['wc','-l','solvate-empty-uncentered.gro'],
				cwd=self.rootdir).split()[0])-3)

		print "translating water box"
		#---removed lipid_water_buffer from below
		lwtranslate = self.settings['solvent_thickness']
		cmd = [gmxpaths['editconf'],
			'-f solvate-empty-uncentered.gro',
			'-o solvate-empty.gro',
			'-translate 0 0 -'+str(lwtranslate)]
		call(cmd,logfile='log-editconf-solvate-move-empty',cwd=self.rootdir)

		print "concatenating water and bilayer"
		#---note that the AAMD setup is more sensitive and likely to crash
		#---...for that reason we use the "editdconf -d 0" output to make it more flush
		#---...and assume that the 
		#---swapped in the md-vacuum-p2-check-dims.gro for vacuum-packed.gro here
		confs = [[line for line in open(self.rootdir+fn,'r')] 
			for fn in [
				('md-vacuum-p2-check-dims.gro'if self.simscale == 'aamd' else 'vacuum-packed.gro'),
				'solvate-empty.gro']]
		natoms = str(sum([len(conf)-3 for conf in confs]))
		fp = open(self.rootdir+'solvate-unsized.gro','w')
		fp.write('bilayer\n'+natoms+'\n')
		for conf in confs:
			for line in conf[2:-1]:
				fp.write(line[:20]+' '+''.join([
					'{value[0]:>3}.{value[1]:<3}'.format(value=str('{0:.2f}'.format(float(i))).split('.'))
					for i in line[20:].strip('\n').split()])+'\n')
		fp.write(confs[0][-1])
		fp.close()

		print "checking the bilayer z-dimension"
		cmd = [gmxpaths['editconf'],
			'-f solvate-unsized.gro',
			'-o solvate-unsized-check-dims.gro',
			'-d '+str(float(self.settings['lipid_water_buffer']))]
		call(cmd,logfile='log-editconf-solvate-unsized',cwd=self.rootdir)
		boxdimsxy = boxdims[:2]
		boxdims = self.get_box_vectors('log-editconf-solvate-unsized',new=True)

		print "resizing the box"
		cmd = [gmxpaths['editconf'],
			'-f solvate-unsized.gro',
			'-o solvate.gro',
			'-box '+str(boxdimsxy[0])+' '+str(boxdimsxy[0])+' '+str(boxdims[2])+' ']
		call(cmd,logfile='log-editconf-solvate',cwd=self.rootdir)

		print "counting waters and updating topology"
		if self.simscale == 'aamd':
			nwaters = (int(checkout(['wc','-l','solvate-empty-uncentered.gro'],
				cwd=self.rootdir).split()[0])-3)/3
		elif self.simscale == 'cgmd':
			nwaters = (int(checkout(['wc','-l','solvate-empty-uncentered.gro'],
				cwd=self.rootdir).split()[0])-3)
		self.lnames.append(self.settings['sol_name'])
		self.comps.append([self.settings['sol_name'],str(nwaters)])
	
		self.write_topology('solvate.top')
		self.minimization_method('solvate')
		
		print "translating so the bilayer is in the middle"
		call('mv solvate-minimized.gro solvate-minimized-unshifted.gro',cwd=self.rootdir)
		cmd = ['echo -e "0\n" |',
			gmxpaths['trjconv'],
			'-f solvate-minimized-unshifted.gro',
			'-o solvate-minimized.gro',
			'-trans 0 0 '+str(-1*self.settings['solvent_thickness']/2.-self.settings['lipid_water_buffer']),
			'-s em-solvate-steep.tpr',
			'-pbc mol']
		call(cmd,logfile='log-trjconv-solvate-shift-center',cwd=self.rootdir)		

	def counterionize(self):
		'''Add counterions to the water slab.'''
		if self.settings['concentration_calc'] == 'exempt_lipids':
			print "running genion on the water only"
			self.write_topology('solvate-water.top',water_only=True)
			cmd = [gmxpaths['grompp'],
				'-f input-em-steep-in.mdp',
				'-po genion-water.mdp',
				'-c solvate-empty.gro',
				'-p solvate-water.top',
				'-o genion-water.tpr']
			call(cmd,logfile='log-grompp-genion-water',cwd=self.rootdir)
			cmd = ['echo -e "keep 0\nr '+self.settings['sol_name']+'\nkeep 1\nq\n" |',
				gmxpaths['make_ndx'],
				'-f solvate-empty.gro',
				'-o counterions-empty.ndx']
			call(cmd,logfile='log-make-ndx-counterions-empty',cwd=self.rootdir)
			print "running diagnostic genion to get ion counts"
			cmd = [gmxpaths['genion'],
				'-s genion-water.tpr',
				'-o counterions-water.gro',
				'-conc '+str(self.settings['ion_strength']),
				'-n counterions-empty.ndx']
			call(cmd,logfile='log-genion-water',cwd=self.rootdir)
			padd1,nadd1 = [int(checkout(["awk","'/Will try to add/ {print $"+str(col)+\
				"}'","log-genion-water"],cwd=self.rootdir).strip()) for col in [5,9]]
			total_ions = padd1+nadd1
			print "checking system charge"
			call('cp solvate.top counterions-charge.top',cwd=self.rootdir)
			cmd = [gmxpaths['grompp'],
				'-f input-em-steep-in.mdp',
				'-po genion-charge.mdp',
				'-c solvate-minimized.gro',
				'-p counterions-charge.top',
				'-o genion-charge.tpr']
			call(cmd,logfile='log-grompp-genion-charge',cwd=self.rootdir)
			cmd = ['echo -e "keep 0\nr '+self.settings['sol_name']+'\nkeep 1\nq\n" |',
				gmxpaths['make_ndx'],
				'-f solvate-minimized.gro',
				'-o counterions-charge.ndx']	
			call(cmd,logfile='log-make-ndx-counterions-charge',cwd=self.rootdir)
			cmd = [gmxpaths['genion'],
				'-s genion-charge.tpr',
				'-o counterions-charge.gro',
				'-nname '+self.settings['salt'][1],
				'-pname '+self.settings['salt'][0],
				'-neutral',
				'-conc '+str(self.settings['ion_strength']),
				'-n counterions-charge.ndx']
			call(cmd,logfile='log-genion-charge',cwd=self.rootdir)
			padd2,nadd2 = [int(checkout(["awk","'/Will try to add/ {print $"+str(col)+\
				"}'","log-genion-charge"],cwd=self.rootdir).strip()) for col in [5,9]]
			charge_difference = padd2-nadd2
			diff = -charge_difference
			total = total_ions
			qn = -self.settings['nq']
			qp = self.settings['pq']
			for t in range(total,total+20):
					p = (qn*t-float(diff))/float(qn+qp)
					n = (qp*t+float(diff))/float(qn+qp)
					if int(p) == float(p) and int(n) == float(n):
						pcount,ncount = int(p),int(n)
						break
			if ncount < 0 or pcount < 0: 
				print "WARNING: zero ions"
			print "adding counterions"
			call('cp solvate.top counterions.top',cwd=self.rootdir)
			self.write_topology('solvate.top')
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
				'-o counterions-waters.ndx']
			call(cmd,logfile='log-make-ndx-counterions-waters',cwd=self.rootdir)
			cmd = [gmxpaths['genion'],
				'-s genion.tpr',
				'-o counterions.gro',
				'-nname '+self.settings['negative_ion_name'],
				'-pname '+self.settings['positive_ion_name'],
				'-nn '+str(ncount),
				'-np '+str(pcount),
				'-n counterions-waters.ndx']
			call(cmd,logfile='log-genion',cwd=self.rootdir)
		if self.settings['concentration_calc'] == 'simple':
			print 'adding counterions'
			call('cp solvate.top counterions.top',cwd=self.rootdir)
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
				'-o counterions-waters.ndx']
			call(cmd,logfile='log-make-ndx-counterions-waters',cwd=self.rootdir)
			cmd = [gmxpaths['genion'],
				'-s genion.tpr',
				'-o counterions.gro',
				'-nname '+self.settings['negative_ion_name'],
				'-pname '+self.settings['positive_ion_name'],
				'-neutral',
				'-conc '+str(self.settings['ion_strength']),
				'-n counterions-waters.ndx']
			call(cmd,logfile='log-genion',cwd=self.rootdir)
		else: raise Exception('except: unclear concentration_calc')

		#---diagnostic
		pcount,ncount = [int(checkout(["awk","'/Will try to add / {print $"+str(col)+\
			"}'","log-genion"],cwd=self.rootdir)) for col in [5,9]]
		if self.simscale == 'aamd':
			nwaters = (int(checkout(['wc','-l','solvate-empty-uncentered.gro'],
				cwd=self.rootdir).split()[0])-3)/3
		elif self.simscale == 'cgmd':
			nwaters = (int(checkout(['wc','-l','solvate-empty-uncentered.gro'],
				cwd=self.rootdir).split()[0])-3)
		self.lnames.append(self.settings['positive_ion_name'])
		self.comps.append([self.settings['positive_ion_name'],str(pcount)])
		self.lnames.append(self.settings['negative_ion_name'])
		self.comps.append([self.settings['negative_ion_name'],str(ncount)])
		self.comps[self.lnames.index(self.settings['sol_name'])][1] = str(nwaters - pcount - ncount)

		self.write_topology('counterions.top')
		self.minimization_method('counterions')

	def grouping(self):
		'''
		Write the ``system-groups.ndx`` file for further simulation, which requires that the water and 
		lipids and possibly proteins are coupled to their own temperature baths.
		'''
		print "writing groups ndx file"
		cmd = ['echo -ne "keep 0\nr '+\
			' | r '.join([l for l in self.lnames if l not in self.ion_residue_names]+['ION'])+\
			'\n'+'r '+' | r '.join([l for l in self.lnames if l in self.ion_residue_names])+\
			'\nname 1 LIPIDS\nname 2 SOLV\ndel 0\nq\n" |',
			gmxpaths['make_ndx'],
			'-f counterions-minimized.gro',
			'-o system-groups.ndx']
		call(cmd,logfile='log-make-ndx-groups',cwd=self.rootdir)
		call('cp counterions-minimized.gro system.gro',cwd=self.rootdir)
		self.write_topology('system.top')

