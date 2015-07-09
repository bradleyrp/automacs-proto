#!/usr/bin/python

#---AUTOMACS (AMX) "LIPIDGRID" MODULE
#---created 2014.08.13 by ryan bradley

"""
This module will generate a lipid monolayer configurations according to input compositions and
a few fairly standard geometric measurements. It relies on a small library of lipids structures and 
topologies which are typically located in the ``sources`` directory.

.. figure:: ../sources/docs/automacs-cgmd-bilayer-monolayer.jpeg
   :align: center
   :alt: monolayer snapshot
   :width: 75 %
   
   Example monolayer configuration.
   
The image above is an example monolayer configuration written to ``place-grid-start.gro``. You can see that 
the coarse-grained beads are not connected. The ``MonolayerGrids`` instance does not generate a topology, 
so we won't see the bonds until we use the ``amx.bilayer`` submodule to construct and minimize the bilayer.
   
Settings
********

Input settings are imported from the python-style file called ``inputs/input_specs_bilayer.py`` by default.
For a bilayer with different compositions for each monolayer, we use the following defintions. ::

	#---COMPOSITIONS
	#---note you can set equal monolayer compositions by deleting complist1, etc

	nlipids0 = 400
	complist0 = [80,19,1]
	lipidnames0 = ['DOPC','DOPS','PIP2']
	nlipids1 = 400
	complist1 = [4,1]
	lipidnames1 = ['DOPC','DOPS']
	
Elsewhere in the ``inputs/input_specs_bilayer.py`` we find the geometric settings for assembling the
monolayers. The key parameter here is the ``solvent_thickness`` which sets the size of the water slab in 
nanometers. ::

	bilayer_construction_settings = {
		'lbuffer':0.1,
		'tbuffer':0.1,
		'pbc_spacing':0.2,
		'solvent_thickness':10,
		'lipid_water_buffer':0.1,
		'ion_strength':0.150,
		'negative_ion_name':'CL','nq':-1,
		'positive_ion_name':'NA','pq':1,
		'lipid_spacing':0.2,
		}
		
These parameters also set the spacing between lipids, between monolayer tails, and so on.

================== ================================================================================
Parameter          Definition
================== ================================================================================
lbuffer            spacing between individual lipids in the grid
tbuffer            spacing between the ends of lipid tails when making a bilayer from a monolayer
pbc_spacing        distance between monolayers under PBC during vacuum construction step
solvent_thickness  thickness of the water slab in nm, normal to the bilayer (may relax)
lipid_water_buffer gap between lipids and water during solvate step
ion_strength       molar ion concentration (check the code for the exact procedure)
lipid_spacing      the distance in nm between adjacent lipids during combinator methods 
================== ================================================================================

"""

#---development features
import os
if os.path.isfile('/etc/pythonstart'):
	execfile('/etc/pythonstart')

#---imports
import os
import sys
import math
import random
from tools import call,checkout,tee

#---locate gmxpaths
if os.path.isfile('gmxpaths.conf'): gmxpaths_path = './'
else: gmxpaths_path = '../'

#---load the gromacs paths
gmxpaths = dict([[i.split()[0],' '.join(i.split()[1:])] 
	for i in open(gmxpaths_path+'gmxpaths.conf','r') 
	if i.strip()[0] != '#'])

#---FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

class MonolayerGrids:
	'''
	Class which generates a lipid monolayer configurations according to desired compositions.
	
	Creating a single instance of this class will automatically run the construction procedure in the root 
	directory (a required keyword argument). Compositions are specified in ``inputs/input_specs_bilayer.py``
	and the user can choose to make the bilayer symmetric by only defining one composition there.

	Parameters
	----------
	rootdir : str
		relative path to perform the monolayer construction and save the final configuration file
	'''
	def __init__(self,rootdir,**kwargs):
		'''
		Sets root directory and passes a dictionary of extra parameters to the build_monolayer function.
		'''
		#---default parameters
		self.params = {
			'structs_dir':'sources/structures/cgmd-bilayer-lipids-structs',
			'input_specs':'inputs/input_specs_bilayer.py',
			}
		#---override defaults with kwargs
		for key in ['structs_dir','input_specs']:
			if key in kwargs.keys(): self.params[key] = kwargs[key]
		self.rootdir = rootdir+'/'

		#---skip everything if the procedure is final configuration is already available
		#---note that this takes the place of an exception and hence expects the user to remove old files
		if os.path.isfile(self.rootdir+'place-grid-start.gro') or \
			os.path.isfile(self.rootdir+'place-grid-start0.gro'):
			sys.stdout = tee(open(self.rootdir+'log-script-master','a',1))
			sys.stderr = tee(open(self.rootdir+'log-script-master','a',1),error=True)
			print 'START LIPIDGRID'
			print 'skipping this step because monolayer configuration exists'
		else:
			#---make root directory
			if not os.path.isdir(os.path.abspath(os.path.expanduser(self.rootdir))): os.mkdir(self.rootdir)
			else: raise Exception('except: requested root directory already exists without a configuration')
			#---start the logger
			sys.stdout = tee(open(self.rootdir+'log-script-master','a',1))
			sys.stderr = tee(open(self.rootdir+'log-script-master','a',1),error=True)
			print 'START LIPIDGRID'
			for key in self.params.keys(): print key+' = '+str(self.params[key])
			#---build the monolayers
			else: self.build_monolayers() 

	def combinator(self,big,add,out,spacing=0.2,rootdir='./'):
		"""
		Add a lipid to a row of lipids.

		This function takes a lipid configuration with an arbitrary number of lipid structures in a straight 
		line. It adds a new lipid to the row, and writes the resulting configuration. It is designed to build 
		a row of lipids piecewise in the y-dimension. The complete rows are assembled by the combinator_row 
		function in the	x-dimension to create a grid of lipids which will eventually become a monolayer. This 
		function is called repeatedly by the :class:`build_monolayers <MonolayerGrids.build_monolayers>` 
		function.
	
		Parameters
		----------
		big : str
			relative path to a structure (gro) file containing a row of lipids
		add : str
			relative path to a structure (gro) file of the lipid to be added to the row of lipids
		out : str
			write the combined configuration to an output structure (gro) file
		"""
	
		call(gmxpaths['editconf']+' -f '+big+' -resnr 2 -o 1.1.gro -center 0 0 0 -d 0',
			logfile='log-editconf-1.1',cwd=rootdir,silent=True)
		call(gmxpaths['editconf']+' -f '+add+' -o 2.1.gro -center 0 0 0 -d 0',
			logfile='log-editconf-2.1',cwd=rootdir,silent=True)
		vy11 = float(checkout(["awk","'/new box vectors/ {print $6}'","log-editconf-1.1"],cwd=rootdir))
		vy21 = float(checkout(["awk","'/new box vectors/ {print $6}'","log-editconf-2.1"],cwd=rootdir))
		call(gmxpaths['editconf']+' -f 1.1.gro -center 0 '+str(vy11/2)+' 0 -o 1.2.gro',
			logfile='log-editconf-1.2',cwd=rootdir,silent=True)
		call(gmxpaths['editconf']+' -f 2.1.gro -center 0 '+str(vy21/2+vy11+spacing)+' 0 -o 2.2.gro',
			logfile='log-editconf-2.2',
			cwd=rootdir,silent=True)
		f = open(rootdir+'3.1.gro', 'w')
		f.write('Test\n')
		n1 = float(checkout(["awk","'NR==2 {print $1}'","1.2.gro"],cwd=rootdir))
		n2 = float(checkout(["awk","'NR==2 {print $1}'","2.2.gro"],cwd=rootdir))
		ntot = n1+n2
		f.write(str(int(ntot))+'\n')
		body1 = ''.join([line for line in open(rootdir+'/1.2.gro')][2:-1])
		body2 = ''.join([line for line in open(rootdir+'/2.2.gro')][2:-1])
		tail = ''.join([line for line in open(rootdir+'/2.2.gro')][-1:])
		f.write(body2)
		f.write(body1)
		f.write(tail)
		f.close()
		call(gmxpaths['editconf']+' -f 3.1.gro -o '+out+' -d 0',
			logfile='log-editconf-out',cwd=rootdir,silent=True)
		call('rm 1.1.gro',cwd=rootdir,silent=True)
		call('rm 2.1.gro',cwd=rootdir,silent=True)
		call('rm 1.2.gro',cwd=rootdir,silent=True)
		call('rm 2.2.gro',cwd=rootdir,silent=True)
		call('rm 3.1.gro',cwd=rootdir,silent=True)
	
	def combinator_row(self,big,add,out,spacing=0.2,rootdir='./'):
		"""
		Add a row of lipids to a grid of lipids.

		This function takes a lipid configuration with an arbitrary number of rows of lipid structures in a 
		grid Formation. It adds a new lipid row to the grid, and writes the resulting configuration. It adds 
		rows in The x-dimension, and hence takes a lipid row which is oriented in the y-dimension. This 
		function is called repeatedly by the build_monolayer function.

		Parameters
		----------
		big : str
			relative path to a structure (gro) file containing a grid of lipids
		add : str
			relative path to a structure (gro) file of the lipid row to be added to the grid of lipids
		out : str
			write the combined configuration to an output structure (gro) file
		"""

		call(gmxpaths['editconf']+' -f '+big+' -o 1.1.gro -center 0 0 0 -d 0',
			logfile='log-editconf-1.1',cwd=rootdir,silent=True)
		call(gmxpaths['editconf']+' -f '+add+' -o 2.1.gro -center 0 0 0 -d 0',
			logfile='log-editconf-2.1',cwd=rootdir,silent=True)
		vy11 = float(checkout(["awk","'/new box vectors/ {print $5}'","log-editconf-1.1"],cwd=rootdir))
		vy21 = float(checkout(["awk","'/new box vectors/ {print $5}'","log-editconf-2.1"],cwd=rootdir))
		call(gmxpaths['editconf']+' -f 1.1.gro -center '+str(vy11/2)+' 0 0 -o 1.2.gro',logfile='log-editconf-1.2',
			cwd=rootdir,silent=True)
		call(gmxpaths['editconf']+' -f 2.1.gro -center '+str(vy21/2+vy11+spacing)+' 0 0 -o 2.2.gro',
			logfile='log-editconf-2.2',cwd=rootdir,silent=True)
		f = open(rootdir+'3.1.gro', 'w')
		f.write('Test\n')
		n1 = float(checkout(["awk","'NR==2 {print $1}'","1.2.gro"],cwd=rootdir))
		n2 = float(checkout(["awk","'NR==2 {print $1}'","2.2.gro"],cwd=rootdir))
		ntot = n1+n2
		f.write(str(int(ntot))+'\n')
		body1 = ''.join([line for line in open(rootdir+'/1.2.gro')][2:-1])
		body2 = ''.join([line for line in open(rootdir+'/2.2.gro')][2:-1])
		tail = ''.join([line for line in open(rootdir+'/2.2.gro')][-1:])
		f.write(body1)
		f.write(body2)
		f.write(tail)
		f.close()
		call(gmxpaths['editconf']+' -f 3.1.gro -o '+out+' -d 0',
			logfile='log-editconf-out',cwd=rootdir,silent=True)
		call('rm 1.1.gro',cwd=rootdir,silent=True)
		call('rm 2.1.gro',cwd=rootdir,silent=True)
		call('rm 1.2.gro',cwd=rootdir,silent=True)
		call('rm 2.2.gro',cwd=rootdir,silent=True)
		call('rm 3.1.gro',cwd=rootdir,silent=True)

	def build_monolayers(self):

		'''
		Generate a random arrangement of lipids and iteratively construct a grid of lipid structures according
		to this arrangement. If ``matplotlib`` is available, the program will render the random configuration
		and save the results in the root directory. Note that this step currently requires ``numpy``.
		'''

		rootdir = self.rootdir
		input_specs = self.params['input_specs']
		lipid_location = os.path.abspath(os.path.expanduser(self.params['structs_dir']))

		#---numpy sqrt below
	
		#---load the compositions
		if os.path.isfile(input_specs):
			#---note that execfile must receive a dictionary because locals cannot be modified (optimization)
			execfile(input_specs,self.params)
			complist0 = self.params['complist0']
			nlipids0 = self.params['nlipids0']
			lipidnames0 = self.params['lipidnames0']
			print 'complist0 = '+str(complist0)
			print 'nlipids0 = '+str(nlipids0)
			print 'lipidnames0 = '+str(lipidnames0)
			totside = [int(math.sqrt(nlipids0))]
			total_lipids = [nlipids0]
		else: raise Exception('except: cannot locate the input specifications file '+str(input_specs))
		#---load compositions for a second, distinct monolayer if available
		if 'complist1' not in self.params.keys(): asymmetric = False
		else:
			print 'detected asymmetric composition from input_specs'
			asymmetric = True
			complist1 = self.params['complist1']
			nlipids1 = self.params['nlipids1']
			lipidnames1 = self.params['lipidnames1']
			print 'complist1 = '+str(complist1)
			print 'nlipids1 = '+str(nlipids1)
			print 'lipidnames1 = '+str(lipidnames1)
			totside = [int(math.sqrt(nlipids0)),int(math.sqrt(nlipids1))]
			total_lipids = [nlipids0,nlipids1]
		self.simscale = self.params['simscale']
		spacing = self.params['bilayer_construction_settings'][self.simscale]['lipid_spacing']
	
		#---record compositions
		if asymmetric: monocomplist,lipidnames = [complist0,complist1],[lipidnames0,lipidnames1]
		else: monocomplist,lipidnames = [complist0],[lipidnames0]
		lipid_nameset = list(set([i for j in lipidnames for i in j]))
		lipid_counts = [0 for i in lipid_nameset]
		
		#---clear previous files to avoid downstream confusion
		for fname in ['./place-grid-start.gro','./place-grid-start0.gro','./place-grid-start1.gro']:
			if os.path.isfile(fname): call('rm '+fname,cwd=rootdir,silent=True)

		#---loop over distinct monolayer compositions	
		for ci in range(len(monocomplist)):
			complist = monocomplist[ci]
			print 'building monolayer '+str(ci)

			#---generate random configuration
			counts = []
			
			if sum(complist) != 1.0: complist = [float(i)/sum(complist) for i in complist]
			arrangement = []
			for i in range(len(complist))[::-1]:
				nlipids = int(round(complist[i]*nlipids0))
				counts.append(nlipids)
				if len(arrangement)+nlipids > total_lipids[ci]: 
					print 'rounding correction on the number of '+lipidnames[ci][i]
					nlipids = total_lipids[ci]-len(arrangement)
				print 'lipid with composition = '+str(complist[i])+' requires '+str(nlipids)+' lipids'
				arrangement.extend([i for x in range(nlipids)])
				lipid_counts[lipid_nameset.index(lipidnames[ci][i])] += nlipids
			
			random.shuffle(arrangement)
			arrangement = [[arrangement[i*totside[ci]+j] for j in range(totside[ci]) 
				if i*totside[ci]+j < len(arrangement)] for i in range(len(arrangement)/totside[ci]+
				(0 if len(arrangement)/totside[ci] == float(len(arrangement))/totside[ci] else 1))]
			print 'proposed configuration = \n'+'\n'.join([str(i) for i in arrangement])

			#---vizualization feature is disabled to prevent errors
			try: 
				import numpy
				import matplotlib as mpl
				import matplotlib.pylab as plt
				cmap = mpl.cm.jet
				fig = plt.figure()
				ax = fig.add_subplot(111)
				ax.set_title('lipid configuration')
				afill = list(arrangement)
				if len(afill[-1]) < len(afill[0]):
					afill[-1] = afill[-1]+[float(len(complist)) 
						for i in range(len(afill[0])-len(afill[-1]))]
					print 'note that fill values of '+str(len(complist))+' are empty'
					bounds = range(len(complist)+1+1)
				else: bounds = range(len(complist)+1)
				afill = numpy.array(afill)
				print bounds
				print afill
				norm = mpl.colors.BoundaryNorm(bounds,cmap.N)
				im = plt.imshow(afill.T,interpolation='nearest',
					origin='lower',cmap=cmap,norm=norm)
				cbar = plt.colorbar(im,orientation="vertical")
				print 'found matplotlib and saving figure to fig-arrangement'+str(ci)+'.png'
				plt.savefig(self.rootdir+'fig-arrangement'+str(ci)+'.png',dpi=100)
				plt.clf()

			except ImportError: print 'no numpy or matplotlib so the text arrangement will have to do'

			#---write the configuration
			fp = open(rootdir+'config.dat','w')
			fp.write('arrangement = '+str(list([list(i) for i in arrangement])))
			fp.close()

			#---generate the grid
			lipidlist = [lipid_location+'/'+i+'.gro' for i in lipidnames[ci]]
			for j in range(0,len(arrangement)):
				start_lipid=''
				next_lipid=''
				start_lipid = lipidlist[arrangement[j][0]]
				next_lipid = lipidlist[arrangement[j][1]]
				self.combinator(start_lipid,next_lipid,'out.gro',spacing=spacing,rootdir=rootdir)	
				for i in range(2,len(arrangement[j])):
					call('mv out.gro in.gro',cwd=rootdir,silent=True)
					next_lipid = lipidlist[arrangement[j][i]]
					self.combinator('in.gro',next_lipid,'out.gro',spacing=spacing,rootdir=rootdir)
				estr=gmxpaths['editconf']+' -f out.gro -o row'+str(j)+'.gro -d 0 -resnr 1'
				call(estr,cwd=rootdir,silent=True,logfile='log-editconf-out')
				call('rm in.gro',cwd=rootdir,silent=True)
				call('rm out.gro',cwd=rootdir,silent=True)
			self.combinator_row('row0.gro','row1.gro','grid.gro',spacing=spacing,rootdir=rootdir)
			for j in range(2,len(arrangement)):
				self.combinator_row('grid.gro','row'+str(j)+'.gro','grid.gro',
					spacing=spacing,rootdir=rootdir)
			if not asymmetric: finalname = 'place-grid-start.gro'
			else: finalname = 'place-grid-start'+str(ci)+'.gro'
			call('cp grid.gro '+finalname,cwd=rootdir,silent=True)
			call('rm row*.gro',cwd=rootdir,silent=True)
			call('rm \#*',cwd=rootdir,silent=True)

		#---write the composition
		fp = open(rootdir+'composition.dat','w')
		for i in range(len(lipid_nameset)): fp.write(lipid_nameset[i]+' '+str(lipid_counts[i])+'\n')	
		fp.close()
		print 'done building monolayer'+('s' if j == 1 else '')
	
