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
import numpy as np
try: from amx.simset import SimSet
except: from simset import SimSet
try:
	from MDAnalysis import *
	import scipy
	from scipy import spatial
	from numpy import where
except: print 'warning: failed to load necessary packages for proteinbilayer.py'

#---locate gmxpaths
if os.path.isfile('gmxpaths.conf'): gmxpaths_path = './'
else: gmxpaths_path = '../'

#---load the gromacs paths
gmxpaths = dict([[i.split()[0],' '.join(i.split()[1:])] 
	for i in open(gmxpaths_path+'gmxpaths.conf','r') 
	if i.strip()[0] != '#'])
	
#---rotation functions
def vecangle(a,b): return np.arccos(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b))
def crossnorm(a,b): return np.cross(a,b)/np.linalg.norm(np.cross(a,b))
def projector(a,dim): return a[np.where(np.arange(3)!=dim)[0]]
def rotation_matrix(axis, theta):
	'''Rotation matrix via Euler-Rodrigues formula.'''
	axis = np.asarray(axis)
	theta = np.asarray(theta)
	axis = axis/np.sqrt(np.dot(axis,axis))
	a = np.cos(theta/2)
	b, c, d = -axis*np.sin(theta/2)
	aa, bb, cc, dd = a*a, b*b, c*c, d*d
	bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
	rotator = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
		[2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
		[2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
	return rotator
	
class ProteinBilayer(amxsim.AMXSimulation):

	def __init__(self,rootdir=None,**kwargs):
		#---set paths for the previous and current steps
		self.rootdir = os.path.abspath(os.path.expanduser(rootdir))+'/'
		
		#---manually specify the sources directory and the inputs file
		if 'sources_dir' in kwargs.keys(): 
			self.sources_dir = os.path.abspath(os.path.expanduser(kwargs['sources_dir']))+'/'
		else: self.sources_dir = os.path.abspath(os.path.expanduser('./sources/'))+'/'
		if 'inputs_file' in kwargs.keys(): 
			self.inputs_file = os.path.abspath(os.path.expanduser(kwargs['inputs_file']))+'/'
		else: 
			self.inputs_file = os.path.abspath(
				os.path.expanduser('./inputs/input-specs-protein-bilayer.dat'))

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
			print 'START PROTEIN + BILAYER'
			#---protein construction settings are pulled from the inputs_file
			self.params = dict()
			execfile(self.inputs_file,self.params)
			self.simscale = self.params['simscale']
			self.settings = self.params['construction_settings'][self.simscale]
			print self.rootdir
			#---copy files from the repo
			copy('sources/cgmd-bilayer-construct/input-em-*-in.mdp',self.rootdir)
			copy('sources/martini.ff',self.rootdir+'martini.ff')
			copy('sources/cgmd-bilayer-lipids-tops',self.rootdir+'lipids-tops')
			for fn in [self.settings[i] for i in 
				['struct_protein','struct_lipid','struct_membrane','top_protein']]:
				copy(fn,self.rootdir+os.path.basename(fn))
			copy(self.settings['top_membrane'],self.rootdir+'system-membrane.top')
			
			
			if 0:
				#---copy input files if this is a fresh start
				copy(self.sources_dir+'aamd-protein-structs/'+self.params['input_filename'],
					self.rootdir+'system-input.pdb')
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
			#---note that this method is only set to simulate a np.single protein
			self.npoz,self.nneg = 0,0
			self.nsol = 0
			self.protname = 'Protein'
			#---call the master construction procedure
			print 'starting protein + bilayer adhesion'
			self.mset = SimSet()
			self.construction()
			self.detect_composition()
			self.remove_close_waters()
			self.convert_ions_to_water()
			copy(self.rootdir+'system-allwater.gro',self.rootdir+'solvate-minimized.gro')
			self.counterionize_general()
			self.minimization_method('counterions')
			call('cp counterions-minimized.gro system.gro',cwd=self.rootdir)
			if not os.path.isfile(self.rootdir+'system.top'): 
				self.write_topology_bilayer('system.top')
			if not os.path.isfile(self.rootdir+'system-groups.ndx'): 
				self.grouping(grouptype='protein-bilayer')
			
	def detect_composition(self):
		
		with open(self.rootdir+'system-membrane.top','r') as fp: 
			top = fp.readlines()
		break1 = [ii for ii,i in enumerate(top) if re.match('^(\s+)?\[(\s+)?molecules(\s+)?\]',i)][0]
		reslist = [[l.split()[0],int(l.split()[1])] for l in top[break1+1:]]
		resnames = [i[0] for i in reslist]
		reslist = [['Protein',self.nprots]]+reslist
		self.lnames,self.comps = [i[0] for i in reslist],reslist
		with open(self.rootdir+'solvate.top','w') as fp:
			break0 = [ii for ii,i in enumerate(top) if re.match('^(\s+)?\[(\s+)?system(\s+)?\]',i)][0]
			break1 = [ii for ii,i in enumerate(top) if re.match('^(\s+)?\[(\s+)?molecules(\s+)?\]',i)][0]
			for li in range(break0-1): fp.write(top[li])
			fp.write('#include "'+self.settings['top_protein'].split('/')[-1]+'"\n')
			fp.write('[ system ]\nPROTEIN AND BILAYER\n[ molecules ]\n')
			reslist = [[l.split()[0],int(l.split()[1])] for l in top[break1+1:]]
			resnames = [i[0] for i in reslist]
			reslist = [['Protein_A',self.nprots]]+reslist
			for r in reslist: fp.write(r[0]+' '+str(r[1])+'\n')
	
	def remove_close_waters(self):
		
		call('editconf -f system-fused.gro -o system-fused-fix.gro -d 0 -resnr 1',
			cwd=self.rootdir,logfile='log-editconf-fix')
		with open(self.rootdir+'script-vmd-remove-waters.tcl','w') as fp:
			fp.write('mol new system-fused-fix.gro'+'\n')
			fp.write('set sel [atomselect top \"all not ((resname W or resname ION) and within ')
			fp.write(str(self.settings['water_gap_angstroms'])+' of (name BB or resname PIP2))\"]\n')
			fp.write('$sel writepdb solvate-vmd.pdb\nexit\n')
		call('vmd -dispdev text -e script-vmd-remove-waters.tcl',
			cwd=self.rootdir,logfile='log-vmd-remove-waters')
		cmd = [gmxpaths['editconf'],
			'-f solvate-vmd.pdb',
			'-o solvate-vmd.gro',
			'-resnr 1']
		call(cmd,logfile='log-editconf-convert-vmd-gro',cwd=self.rootdir)
		
	def convert_ions_to_water(self):
		
		with open(self.rootdir+'solvate-vmd.gro','r') as fp: gro = fp.readlines()
		with open(self.rootdir+'system-allwater.gro','w') as fp:
			fp.write(gro[0])
			fp.write(gro[1])
			for line in gro[2:-1]:
				for iname in ['negative_ion_name','positive_ion_name']:
					line = re.sub('ION'+(re.sub('\+','\\+',
						self.settings[iname]).rjust(7+(1 if '+' in self.settings[iname] else 0))),
						'W'.ljust(5)+'W'.rjust(5),line)
				fp.write(line)
			fp.write(gro[-1])

		print "counting waters and updating topology"
		with open(self.rootdir+'system-allwater.gro','r') as fp: gro = fp.readlines()
		nwaters = len([i for i in gro if re.search('W'.ljust(5)+'W'.rjust(5),i)])
		print self.lnames,self.comps
		delinds = [self.lnames.index(self.settings[i]) 
			for i in ['sol_name','negative_ion_name','positive_ion_name']]
		print np.sort(delinds)
		for i in np.sort(delinds)[::-1]: 
			del self.lnames[i]
			del self.comps[i]
		self.lnames.append('W')
		self.comps.append(['W',nwaters])
		print self.lnames,self.comps
		self.protein_itp = os.path.basename(self.settings['top_protein'])
		self.write_topology_bilayer('solvate.top')

	def construction(self):
	
		#---input specs
		basedir = self.rootdir	
		struct_membrane = os.path.basename(self.settings['struct_membrane'])
		director_cgmd = self.settings['director_cgmd']
		selector_cgmd = self.settings['selector_cgmd']
		top_membrane = self.settings['top_membrane']
		top_protein = self.settings['top_protein']
		struct_protein = os.path.basename(self.settings['struct_protein'])
		method_lattice = self.settings['method_lattice']
		noligos_square = self.settings['noligos_square']
		noligos_hex = self.settings['noligos_hex']
		spacer = self.settings['spacer']
		jitter = self.settings['jitter']
		zoffset = self.settings['zoffset']
		protein_key_residues = self.settings['protein_binding_pocket_residues']

		#---read membrane topology to get a residue list
		#---note that this is a generic topology reader
		fn = top_membrane
		with open(fn,'r') as fp: lines = fp.readlines()
		toplist = [[k[0],int(k[1])] for k in [l.strip('\n').split() 
			for l in lines[[ii for ii,i in enumerate(lines) 
			if re.match('^(\s+)?\[(\s+)?molecules',i)][0]+1:]]]
		solvent_names = [self.settings[i] for i in ['negative_ion_name','positive_ion_name','sol_name']]
		residues = [i for i in [k[0] for k in toplist] if i not in solvent_names]
		
		self.mset.load_trajectory((basedir+struct_membrane,basedir+struct_membrane),resolution='cgmd')
		self.mset.identify_monolayers(director_cgmd,startframeno=0)
		self.mset.identify_residues(residues)
		topxyz = np.array([np.mean(
			self.mset.universe.residues[i].selectAtoms(selector_cgmd).coordinates(),axis=0) 
			for i in self.mset.monolayer_residues[0]])
		vecs = self.mset.vec(0)

		#---load the bilayer gro file
		with open(basedir+struct_membrane,'r') as fp: grobilayer = fp.readlines()

		#---load the protein gro file and read positions
		with open(basedir+struct_protein,'r') as fp: grolines = fp.readlines()
		groprotein = [l for l in grolines[2:-1] if l[10:15].strip() not in solvent_names]
		protpos = np.array([[float(j) for j in i[20:].strip('\n').split()] for i in groprotein])

		#---formatting trick to ensure aligned decimal places
		dotplace = lambda n: re.compile(r'(\d)0+$').sub(r'\1',"%8.2f"%float(n))

		#---method in which we specify spacing between the oligomer and put the lattice 
		#---...in the center of the box formerly called 'inter-oligo-spacing' but now it is the default
		#---generate lateral lattice positions
		#---place oligomers in a square
		if method_lattice == 'square':
			center_lattice = [[i*spacer,j*spacer] 
				for j in range(noligos_square[1]) for i in range(noligos_square[0])]
			corner = vecs[:2]/2. - 10*np.mean(center_lattice,axis=0) + jitter
			pos_oligos = corner + 10*np.array(center_lattice) 
		#---place oligomers in a hexagonal lattice according to the provided list
		if method_lattice == 'hex':
			nx,ny = max(noligos_hex),len(noligos_hex)
			big_dims = 10.*np.array([spacer*(nx-1),np.sqrt(2)/2.*spacer*(ny-1)])
			maxcols = [max([noligos_hex[i] for i in range(k,len(noligos_hex),2)]) for k in range(2)]
			corner = vecs[:2]/2. - big_dims/2. -\
				([spacer*10/4.,0.] if maxcols[0] == maxcols[1] else [0.,0.]) + jitter
			pos_oligos = [[10*spacer*ind1+corner[0]+\
				(10*spacer/2. if noligos_hex[ind2] <= max(noligos_hex) and ind2%2 == 1 else 0.)+\
				(-10*spacer/2. if noligos_hex[ind2] > max(noligos_hex) and ind2%2 == 1 else 0.),
				10*spacer*ind2*np.sqrt(2)/2.+corner[1]] 
				for ind2 in range(len(noligos_hex)) for ind1 in range(noligos_hex[ind2])]
		#---generate monomer positions within each oligomer
		#---note that we previously had an option for nmonos==1 here but it is now only the default
		pos_lateral = pos_oligos

		#---find appropriate heights
		cd = scipy.spatial.distance.cdist(pos_lateral,topxyz[:,:2])
		pos_norms = [np.mean(topxyz[([where(cd[j]==i)[0][0] 
			for i in np.sort(cd[j])[:10]])][:,2])	for j in range(len(cd))]

		#---determine position of the input molecules in order to perform the shifts
		protpos_com = np.mean(protpos[np.array(protein_key_residues)],axis=0)
		protpos_com_true = np.mean(protpos,axis=0)
		centers = np.concatenate((pos_oligos.T,np.array([pos_norms])+10*zoffset)).T
		natoms = (len(groprotein)-3)*len(centers)+len(grobilayer)-3
		self.nprots = len(centers)
		
		#---rotate the protein
		resids = np.array([int(i[:5]) for i in grolines[2:-1]])
		resnames = np.array([i[5:10].strip() for i in grolines[2:-1]])
		atomnames = np.array([i[10:15].strip() for i in grolines[2:-1]])
		keypos = protpos[np.array([np.where(np.all(np.array([resids==i,
			atomnames=='BB']).T,axis=1))[0][0] for i in protein_key_residues])]
		keycom = np.array([np.mean(keypos,axis=0)])
		allcom = np.array([np.mean(protpos,axis=0)])
		if self.settings['method_adhesion'] == 'ready':
			protpos_final = protpos-np.mean(keypos,axis=0)
		elif self.settings['method_adhesion'] == 'contact':
			rotmat = rotation_matrix(crossnorm((keycom-allcom)[0],np.array([0,0,-1])),
				1*vecangle((keycom-allcom)[0],np.array([0,0,-1])))
			protpos_rotated = np.dot((protpos-np.mean(protpos,axis=0)),rotmat.T)
			#---re-center the key atom center of mass for translating
			keypos2 = np.mean(protpos_rotated[np.array([np.where(np.all(np.array([resids==i,
				atomnames=='BB']).T,axis=1))[0][0] for i in protein_key_residues])],axis=0)
			protpos_final = protpos_rotated-keypos2
		
		#---lipid to delete
		del_list = []
		
		#---write the proteins
		fpout = open(basedir+'/system-fused.gro','w')
		fpout.write('protein+bilayer\n')
		fpout.write(str(natoms)+'\n')
		#---for each protein location print the rotated, shifted coordinates
		for ind in range(len(centers)):
			for gi,line in enumerate(groprotein):
				fpout.write(line[:20]+''.join([dotplace(i) 
					for i in protpos_final[gi]+centers[ind]/10.])+'\n')
		#---write the bilayer after the proteins
		for l in range(2,len(grobilayer)-1):
			if any([(True if re.search(str(self.mset.monolayer_residues[0][del_inds[i]]).rjust(5)+\
				replace_resname,grobilayer[l]) != None else False) 
				for i in del_list]) == False:
				fpout.write(grobilayer[l])
		fpout.write(grobilayer[-1])
		fpout.close()
			
