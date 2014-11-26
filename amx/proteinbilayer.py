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
from amx.simset import SimSet
import amxsim
import numpy as np
from MDAnalysis import *
import scipy
from scipy import spatial
from numpy import where

#---locate gmxpaths
if os.path.isfile('gmxpaths.conf'): gmxpaths_path = './'
else: gmxpaths_path = '../'

#---load the gromacs paths
gmxpaths = dict([[i.split()[0],' '.join(i.split()[1:])] 
	for i in open(gmxpaths_path+'gmxpaths.conf','r') 
	if i.strip()[0] != '#'])
	
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
			#---copy files from the repo
			copy('./repo/*',self.rootdir)
			copy('./sources/cgmd-bilayer-construct/input-em-*-in.mdp',self.rootdir)
			
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
			self.nprots = 1
			self.npoz,self.nneg = 0,0
			self.nsol = 0
			self.protname = 'Protein'
			#---call the master construction procedure
			print 'starting protein + bilayer adhesion'
			self.mset = SimSet()
			self.construction()
			self.minimization_method('counterions')
			
	def construction(self):
	
		'''
		preparation steps

		0. prepare only a lipid and protein gro file, it doesn't matter where they are absolutely located

		steps to create protein lattice

		1. load protein, lipid, and membrane
		2. define the lattice relative to the box vectors
		3. assign lateral positions and rotations
		4. find nearest several proteins and take their average height
		5. specify the height of the complex
		6. specify the offset and determine the z-shift
		7. print or write the x, y, z, and rotation coordinates to the screen or somehow send to bash
		8. use bash or call bash functions to actually write the gro files
		9. once the complex is assembled, start the membed procedure
		10. consider outputting the index numbers for membed

		remaining tasks
		find the nearest DOPS molecule for each PIP2 and delete it before writing the file
		interface directly with g_membed
		'''

		#---used re module to write this code to unpack the dictionary because locals intransigent
		director_cgmd = self.settings['director_cgmd']
		selector_cgmd = self.settings['selector_cgmd']
		struct_protein = self.settings['struct_protein']
		struct_lipid = self.settings['struct_lipid']
		struct_membrane = self.settings['struct_membrane']
		top_membrane = self.settings['top_membrane']
		top_protein = self.settings['top_protein']
		method = self.settings['method']
		method_lattice = self.settings['method_lattice']
		method_lattice = self.settings['method_lattice']
		spacer = self.settings['spacer']
		noligos_square = self.settings['noligos_square']
		noligos_hex = self.settings['noligos_hex']
		spacer_within = self.settings['spacer_within']
		natoms_replace = self.settings['natoms_replace']
		nmonos = self.settings['nmonos']
		rotation_list = self.settings['rotation_list']
		protein_key_residues = self.settings['protein_key_residues']
		replace_resname = self.settings['replace_resname']
		residues = self.settings['residues']
		jitter = self.settings['jitter']
		zoffset = self.settings['zoffset']
		basedir = self.rootdir
		
		self.mset.load_trajectory((basedir+struct_membrane,basedir+struct_membrane),resolution='cgmd')
		self.mset.identify_monolayers(director_cgmd,startframeno=0)
		self.mset.identify_residues(residues)
		topxyz = np.array([np.mean(
			self.mset.universe.residues[i].selectAtoms(selector_cgmd).coordinates(),axis=0) 
			for i in self.mset.monolayer_residues[0]])
		vecs = self.mset.vec(0)

		#---load the protein
		prot = Universe(basedir+struct_protein,basedir+struct_protein)
		protpos = prot.selectAtoms('all').coordinates()

		#---load the lipid
		lipid = Universe(basedir+struct_lipid,basedir+struct_lipid)
		lipidpos = lipid.selectAtoms('all').coordinates()

		#---load the bilayer gro file
		fpin_bilayer = open(basedir+struct_membrane,'r')
		struct_bilayer_file = []
		for line in fpin_bilayer:
			struct_bilayer_file.append(line)
		fpin_bilayer.close()

		#---load the protein gro file
		fpin = open(basedir+struct_protein,'r')
		struct_protein_file = []
		for line in fpin:
			struct_protein_file.append(line)
		fpin.close()

		#---load the protein gro file
		fpin = open(basedir+struct_lipid,'r')
		struct_lipid_file = []
		for line in fpin:
			struct_lipid_file.append(line)
		fpin.close()

		#---method in which we specify spacing between the oligomer and put the lattice 
		#---...in the center of the box
		if method == 'inter-oligo-spacing':
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
					for ind2 in range(len(noligos_hex))for ind1 in range(noligos_hex[ind2])]
			#---generate monomer positions within each oligomer
			if nmonos == 1:
				pos_lateral = pos_oligos
				#---find appropriate heights
				cd = scipy.spatial.distance.cdist(pos_lateral,topxyz[:,:2])
				pos_norms = [np.mean(topxyz[([where(cd[j]==i)[0][0] 
					for i in np.sort(cd[j])[:10]])][:,2])	for j in range(len(cd))]
				#---determine position of the input molecules in order to perform the shifts
				protpos_com = np.mean(protpos[np.array(protein_key_residues)],axis=0)
				lipidpos_com = np.mean(lipidpos,axis=0)
				protpos_com_true = np.mean(protpos,axis=0)
				#---write the combined gro file
				fpout = open(basedir+'/system-fused.gro','w')
				fpout.write('protein\n')
				fpout.write(str(len(pos_lateral)*(len(struct_protein_file)-3)+\
					len(pos_lateral)*(len(struct_lipid_file)-3)+\
					len(struct_bilayer_file)-3-\
					natoms_replace*len(pos_lateral))+'\n')
				#---loop over monomers
				del_list = []
				for ind in range(len(pos_lateral)):
					#---compute lateral shift
					shifter = np.concatenate((pos_lateral[ind],[pos_norms[ind]])) - protpos_com
					#---perform rotations along the z axis if desired
					if rotation_list != None:
						theta = rotation_list[ind]
						rotate_mat = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
						rotated = np.array([np.concatenate((rotate_mat.dot(i[:2]),[i[2]])) 
							for i in (protpos-protpos_com_true)])
						out_coords = rotated + protpos_com_true + shifter
					else: out_coords = protpos + shifter
					#---output protein points
					for l in range(2,len(struct_protein_file)-1):
						fpout.write(struct_protein_file[l][0:22]+''+\
							str('%.3f'%(out_coords[l-2][0]/10.))+'  '+\
							str('%.3f'%(out_coords[l-2][1]/10.))+'  '+\
							str('%.3f'%(out_coords[l-2][2]/10.+\
							zoffset))+'\n')
				for ind in range(len(pos_lateral)):
					#---output shifted accompanying lipids
					shifter_lipid = np.concatenate((pos_lateral[ind],[pos_norms[ind]])) +\
						lipidpos_com - protpos_com
					#---perform rotations along the z axis if desired
					if rotation_list != None:
						theta = rotation_list[ind]
						rotate_mat = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
						rotated = np.array([np.concatenate((rotate_mat.dot(i[:2]),[i[2]])) 
							for i in (lipidpos-lipidpos_com)])
						out_coords = rotated + lipidpos_com + shifter_lipid
					else: out_coords = lipidpos + shifter_lipid
					for l in range(2,len(struct_lipid_file)-1):
						fpout.write(struct_lipid_file[l][0:22]+''+str('%.3f'%(out_coords[l-2][0]/10.))+'  '+
							str('%.3f'%(out_coords[l-2][1]/10.))+'  '+str('%.3f'%(out_coords[l-2][2]/10.+\
							zoffset))+'\n')
					#---find nearest lipid to delete and add it to the delete list
					del_inds = [self.mset.monolayer_residues[0].index(i) 
						for i in self.mset.monolayer_by_resid_abs[0][
							self.mset.resnames.index(replace_resname)]]
					del_list.append(np.array([np.linalg.norm(np.mean(out_coords,axis=0)-i) 
						for i in topxyz[del_inds]]).argmin())
				#---write the bilayer
				for l in range(2,len(struct_bilayer_file)-1):
					if any([(True if re.search(str(self.mset.monolayer_residues[0][del_inds[i]]).rjust(5)+\
						replace_resname,struct_bilayer_file[l]) != None else False) 
						for i in del_list]) == False:
						fpout.write(struct_bilayer_file[l])
				fpout.write(struct_bilayer_file[-1])
				fpout.close()


	
