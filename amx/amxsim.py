#!/usr/bin/python

"""
This module contains the parent class which is generalized across almost all GROMACS procedures in automacs.
"""

from tools import call,checkout,tee,copy,chain_steps,latestcheck,lastframe,bp,delve,mapdict,redundant
import os,sys,re

#---locate gmxpaths
if os.path.isfile('gmxpaths.conf'): gmxpaths_path = './'
else: gmxpaths_path = '../'

#---load the gromacs paths
gmxpaths = dict([[i.split()[0],' '.join(i.split()[1:])] 
	for i in open(gmxpaths_path+'gmxpaths.conf','r') 
	if i.strip()[0] != '#'])
	
#---CLASS
#-------------------------------------------------------------------------------------------------------------

class AMXSimulation:

	"""
	Parent class for many different types of GROMACS simulations. 
	Comprised of common functions necessary for building most systems.
	"""

	def get_box_vectors(self,logfile,new=True):

		"""
		Read box vectors from editconf output stored in a log file.
		"""

		if new:
			boxdims = [float(checkout(["awk","'/new box vectors/ {print $"+\
				str(i+5)+"}'",logfile],cwd=self.rootdir).strip()) for i in range(3)]
		else:
			boxdims = [float(checkout(["awk","'/    box vectors/ {print $"+\
				str(i+4)+"}'",logfile],cwd=self.rootdir).strip()) for i in range(3)]
		return boxdims
		
	def write_topology_bilayer(self,topname,water_only=False):

		"""
		Write the bilayer topology file from residue lists stored in ``self.lnames`` and ``self.complist``.
		"""

		print 'writing topology'
		#---atomistic topology disabled here
		if self.simscale == 'aamd':
			fp = open(self.rootdir+topname,'w')
			for itp in self.settings['top_includes']: fp.write('#include "'+itp+'"\n')
			if hasattr(self,'protein_itp'): fp.write('#include "'+self.protein_itp+'"\n')
			for lname in self.lnames:
				if lname not in self.ion_residue_names: 
					fp.write('#include "./lipids-tops/lipid.'+lname+'.itp"\n')
			fp.write('[ system ]\nBILAYER\n\n[ molecules ]\n')
			for l in range(len(self.lnames)): 
				if not water_only or self.lnames[l] == self.settings['sol_name']:
					fp.write(self.lnames[l]+' '+self.comps[l][1]+'\n')
			fp.close()
		elif self.simscale == 'cgmd':
			#---define ion names for exclusion from topology update function
			self.ion_residue_names = [self.settings['sol_name'],'MG','NA','CL','Cal'+'ION'+'Ion']+\
				[self.settings[i] for i in ['negative_ion_name','positive_ion_name']]
			fp = open(self.rootdir+topname,'w')
			for itp in self.settings['top_includes']: fp.write('#include "'+itp+'"\n')
			if hasattr(self,'protein_itp'): fp.write('#include "'+self.protein_itp+'"\n')
			fp.write('[ system ]\nBILAYER\n\n[ molecules ]\n')
			print self.comps
			print self.lnames
			for l in range(len(self.lnames)): 
				if not water_only or self.lnames[l] == self.settings['sol_name']: 
					fp.write(self.lnames[l]+' '+str(self.comps[l][1])+'\n')
			fp.close()
		else: raise Exception('except: unclear procedure to write topology')

	def write_topology_protein(self,topname):

		"""
		Write a protein+water topology file.
		"""
		
		print "printing protein topology file"
		if self.simscale == 'aamd':
			fp = open(self.rootdir+topname,'w')
			for itp in self.settings['top_includes']: fp.write('#include "'+itp+'"\n')
			if hasattr(self,'protein_itp'): fp.write('#include "'+self.protein_itp+'"\n')
			fp.write('[ system ]\nprotein+water\n\n[ molecules ]\n')
			if type(self.protname) == list:
				for i in range(len(self.protname)):
					fp.write(self.protname[i]+' '+str(self.nprots[i])+'\n')
			else: fp.write(self.protname+' '+str(self.nprots)+'\n')
			fp.write(self.settings['sol_name']+' '+str(self.nsol)+'\n')
			fp.write(self.settings['positive_counterion_name']+' '+str(self.npoz)+'\n')
			fp.write(self.settings['negative_counterion_name']+' '+str(self.nneg)+'\n')
			fp.close()
		if self.simscale == 'cgmd':
			fp = open(self.rootdir+topname,'w')
			for itp in self.settings['top_includes']: fp.write('#include "'+itp+'"\n')
			if hasattr(self,'protein_itp'): fp.write('#include "'+self.protein_itp+'"\n')
			fp.write('[ system ]\nBILAYER\n\n[ molecules ]\n')
			for lipid_itp in self.itp_lipid: fp.write('#include "lipids-tops/'+lipid_itp+'.itp"'+'\n')
			for protein_itp in self.itp_protein: fp.write('#include "'+protein_itp+'.itp"'+'\n')
			fp.write('[ system ]'+'\n')
			fp.write('PROTEIN, MARTINI, IN WATER'+'\n')
			fp.write('[ molecules ]'+'\n')
			for i in range(len(self.itp_protein)): 
				fp.write(self.itp_protein[i]+' '+str(self.nprots[i])+'\n')
			fp.write(self.settings['sol_name']+' '+str(self.nsol)+'\n')
			fp.write(self.settings['positive_counterion_name']+' '+str(self.npoz)+'\n')
			fp.write(self.settings['negative_counterion_name']+' '+str(self.nneg)+'\n')
			fp.close()
	
	def minimization_method(self,name,posre=False):

		"""
		Generic minimization method with a drop-in name and standard inputs. This function takes a file
		suffix and performs a minimization on the corresponding ``gro`` file.
		"""
		
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
		
	def grouping(self,grouptype='standard',startstruct='counterions-minimized.gro'):

		"""
		Write the ``system-groups.ndx`` file for further simulation, which requires that the water and 
		lipids and possibly proteins are coupled to their own temperature baths.
		"""

		print "writing groups ndx file"
		if grouptype == 'standard':
			#---write a naive groups file
			cmd = [gmxpaths['make_ndx'],
				'-f system.gro',
				'-o system-groups.ndx']
			call(cmd,logfile='log-make-ndx-groups',cwd=self.rootdir,inpipe='q\n')
		elif grouptype == 'cgmd_water':
			inpipe = 'keep 1\nr W | r ION\nname 0 PROTEIN\nname 1 SOLV\nq\n'
			#---write a groups file suitable for bilayers with separate lipid and solvent groups
			cmd = [gmxpaths['make_ndx'],
				'-f counterions-minimized.gro',
				'-o system-groups.ndx']
			call(cmd,logfile='log-make-ndx-groups',cwd=self.rootdir,inpipe=inpipe)
		elif grouptype == 'bilayer':
			inpipe = 'keep 0\nr '+\
				' | r '.join([l for l in self.lnames if l not in self.ion_residue_names+['ION']])+\
				'\n'+'r '+' | r '.join([l for l in self.lnames if l in self.ion_residue_names]+['ION'])+\
				'\nname 1 LIPIDS\nname 2 SOLV\ndel 0\nq\n'
			#---write a groups file suitable for bilayers with separate lipid and solvent groups
			cmd = [gmxpaths['make_ndx'],
				'-f '+startstruct,
				'-o system-groups.ndx']
			call(cmd,logfile='log-make-ndx-groups',cwd=self.rootdir,inpipe=inpipe)
		elif grouptype == 'protein-bilayer':
			inpipe = 'keep 1\nr '+\
				' | r '.join([l for l in self.lnames if l not in self.ion_residue_names+['ION'] and
				not re.match('(PROTEIN|Protein)',l)])+\
				'\n'+'r '+' | r '.join([l for l in self.lnames if l in self.ion_residue_names]+['ION'])+\
				'\nname 0 PROTEIN\nname 1 LIPIDS\nname 2 SOLV\nq\n'
			#---write a groups file suitable for bilayers with separate lipid and solvent groups
			cmd = [gmxpaths['make_ndx'],
				'-f '+startstruct,
				'-o system-groups.ndx']
			call(cmd,logfile='log-make-ndx-groups',cwd=self.rootdir,inpipe=inpipe)
		else: raise Exception('except: incomprehensible group type')
		
	def counterionize_general(self):

		"""
		Add counterions to the water slab.
		"""

		print 'adding counterions'
		call('cp solvate.top counterions.top',cwd=self.rootdir)
		cmd = [gmxpaths['grompp'],
			'-f input-em-steep-in.mdp',
			'-po genion.mdp',
			'-c solvate-minimized.gro',
			'-p counterions.top',
			'-o genion.tpr']
		call(cmd,logfile='log-grompp-genion',cwd=self.rootdir)
		cmd = [gmxpaths['make_ndx'],
			'-f solvate-minimized.gro',
			'-o counterions-waters.ndx']
		call(cmd,logfile='log-make-ndx-counterions-waters',cwd=self.rootdir,
			inpipe='keep 0\nr '+self.settings['sol_name']+'\nkeep 1\nq\n')
		cmd = [gmxpaths['genion'],
			'-s genion.tpr',
			'-o counterions.gro',
			'-nname '+self.settings['negative_ion_name'],
			'-pname '+self.settings['positive_ion_name'],
			'-neutral',
			'-conc '+str(self.settings['ion_strength']),
			'-n counterions-waters.ndx']
		call(cmd,logfile='log-genion',cwd=self.rootdir)
	
		#---diagnostic
		pcount,ncount = [int(checkout(["awk","'/Will try to add / {print $"+str(col)+\
			"}'","log-genion"],cwd=self.rootdir)) for col in [5,9]]
		#---note the following is cgmd only since this routine is set up for proteinbilayer.py
		with open(self.rootdir+'solvate-minimized.gro','r') as fp: gro = fp.readlines()
		nwaters = len([i for i in gro if re.search('W'.ljust(5)+'W'.rjust(5),i)])
		self.lnames.append(self.settings['positive_ion_name'])
		self.comps.append([self.settings['positive_ion_name'],str(pcount)])
		self.lnames.append(self.settings['negative_ion_name'])
		self.comps.append([self.settings['negative_ion_name'],str(ncount)])
		self.comps[self.lnames.index(self.settings['sol_name'])][1] = str(nwaters - pcount - ncount)

		self.write_topology_bilayer('counterions.top')
		self.minimization_method('counterions')

	def resituate(self,basename,tpr):

		"""
		After minimization steps it is often useful to stop lipid from jumping across the box.
		"""
		
		copy(self.rootdir+basename+'.gro',self.rootdir+basename+'-jumped.gro')
		call('rm '+basename+'.gro',cwd=self.rootdir)
		cmd = [gmxpaths['trjconv'],
			'-f '+basename+'-jumped.gro',
			'-o '+basename+'.gro',
			'-s '+tpr+'.tpr',
			'-pbc nojump']
		call(cmd,logfile='log-trjconv-nojump-'+basename,cwd=self.rootdir,inpipe="0\n")				
		
	def write_mdp(self,mdp_section='mdp',rootdir=None):
	
		"""
		Universal MDP file writer which creates input files based on a unified dictionary.
		All MDP parameters should be stored in *mdpdict* within ``inputs/input-specs-mdp.dat``. 
		We assemble MDP files according to rules that can be found in the *mdp* entry for the specs 
		dictionary contained in the ``input-specs-PROCEDURE.dat`` file.
		
		The *mdpdict* variable is a nested dictionary which contains individual chunks of an MDP file. 
		This function parses *mdpdict* and writes the terminal dictionaries (the leaves in a large tree) in
		the standard GROMACS format according to a set of rules which are designed to make it easy to specify
		a sequence of parameters with minimal redundancy. We describe the rules according to the *route*,
		a list of dictionary keys, necessary to navigate from the top of *mdpdict* to a terminal dictionary 
		(a leaf).
		
		1. The *mdp* entry in ``input-specs-PROCEDURE.dat`` contains three types of entries.

			a. The ``top`` entry whittles *mdpdict* according to a list of keywords.
			b. The ``defaults`` entry specifies the route to parameters which should be included in all MDP files.
			c. Any entries that contain the name of an MDP file should point to a list of tuples, where each tuple provides a list of routes necessary to construct the file according to the remaining rules in this list.

		2. A tuple which contains an incomplete route (that is, one that ends in a sub-tree instead of a leaf forces the writer to include all of the entries in that sub-tree and exclude other entries with routes that end in the same key.
		3. If there are any routes that are identical except for the last key, the user must specify the route to the option that they prefer, otherwise an exception will be raised.
		4. Otherwise, the entire (whittled) *mdpdict* is included in the MDP file. This allows the writer to infer default values. These default values can be overridden by placing a leaf with the same name at a lower level of the tree, specifically inside of a sub-tree.
		
		The author welcomes any advice to make these rules more clear.
		"""
	
		target_directory = self.rootdir if rootdir == None else rootdir
		mdpspecs = self.params['bilayer_construction_settings']['cgmd'][mdp_section]
		mdpfile = {}
		execfile('./inputs/input-specs-mdp.dat',mdpfile)
		mdpdict = mdpfile['mdpdict']		

		#--topkeys is the root node for our parameters in mdpdict
		topkeys = mdpspecs['top'] if 'top' in mdpspecs else ()
		#---defspecs set any default paths on the mdpdict tree which we use for all mdp files
		defspecs = mdpspecs['defaults'] if 'defaults' in mdpspecs else ()
		for mdpname in [i for i in mdpspecs if re.match('.+\.mdp$',i)]:
			print 'generating mdp: '+mdpname
			#---for each mdp file the specs defined in the dictionary help select between redundant options
			specs = mdpspecs[mdpname]+defspecs
			routes = list(mapdict(delve(mdpdict,*topkeys)))
			#---divide specs into fulls and partials
			partials = [spec for spec in specs if spec not in routes]
			fulls = [spec for spec in specs if spec in routes]
			#---for each partial we get all of the downstream dictionaries
			downs = [r for r in routes if any([r[:len(p)]==p for p in partials])]
			#---whittle all routes so that they don't include any keyword overlap with downs
			distincts = [r for r in routes if not any([any([i in r for i in d]) for d in downs])]
			#---resolve conflicts with any fully-specified routes
			resolve = [r for r in distincts+downs if not any([r[:-1]==f[:-1] for f in fulls])]
			possibles = resolve+fulls
			#---finally we override all short routes with longer ones
			finals = list(set([v[-1] for v in possibles]))
			check_redundant = [sum([w==f for w in [v[-1] for v in possibles]]) for f in finals]
			dups = [ii for ii,i in enumerate(check_redundant) if i>1]
			valids = [a[0] if len(a)==1 else 
				[b for b in a if len(b)==max([len(c) for c in a])][0] 
				for a in [[v for v in possibles if v[-1]==f] 
				for f in finals]]
			keycollect = [delve(mdpdict,*(topkeys+r)).keys() for r in valids]
			if redundant([i for j in keycollect for i in j]): 
				raise Exception('repeated mdp keys for '+mdpname+': '+str(valids))
			with open(target_directory+'/'+mdpname,'w') as fp:
				for v in valids:
					for key,val in delve(mdpdict,*(topkeys+v)).items():
						fp.write(str(key)+' = '+str(val)+'\n')

		if 0: 
			#--topkeys is the root node for our parameters in mdpdict
			topkeys = mdpspecs['top'] if 'top' in mdpspecs else ()
			#---defspecs set any default paths on the mdpdict tree which we use for all mdp files
			defspecs = mdpspecs['defaults'] if 'defaults' in mdpspecs else ()
			for mdpname in [i for i in mdpspecs if re.match('.+\.mdp$',i)]:
				print 'generating mdp: '+mdpname
				#---for each mdp file the specs defined in the dictionary help select between redundant options
				specs = mdpspecs[mdpname]+defspecs
				routes = list(mapdict(delve(mdpdict,*topkeys)))
				#---any incomplete routes (partials) imply inclusion of all routes that match
				partials = [spec for spec in specs if spec not in routes]
				fulls = [spec for spec in specs if spec in routes]
				"""
				rules for selecting routes
					1. Include all full routes
					2. If a route ends before the terminus of the mdpdict tree keep all the rules below it.
					3. Then discard any rules that are identical to rules in specs except in the final item.
					4. Include any other terminal dictionaries in mdpdict as long as they have no overlapping 
						keys with the rules that lie beneath the partial rules.
				This minimal rule set means that you can exclude large portions of the dictionary by providing a 
				partial rule and making sure that you use common keywords in the sub-rules so that they will 
				exclude each other. This also means that you can override defaults by burying new 
				sub-dictionaries inside of these partial paths.
				"""
				valids_partials = [r for r in routes if 
					(any([r[:len(p)]==p for p in partials]) and
					r[:-1] not in [s[:-1] for s in specs])
					]
				valids = valids_partials + fulls + [r for r in routes if 
					not any([r[:-1]==f[:-1] for f in fulls]) and
					not any([any([v in r for v in vp]) for vp in valids_partials])
					]
				keycollect = [delve(mdpdict,*(topkeys+r)).keys() for r in valids]
				if redundant([i for j in keycollect for i in j]): 
					raise Exception('repeated mdp keys: '+str(valids))
				with open(mdpname,'w') as fp:
					for v in valids:
						for key,val in delve(mdpdict,*(topkeys+v)).items():
							fp.write(str(key)+' = '+str(val)+'\n')

		if 0:
			#--topkeys is the root node for our parameters in mdpdict
			topkeys = mdpspecs['top'] if 'top' in mdpspecs else ()
			#---defspecs set any default paths on the mdpdict tree which we use for all mdp files
			defspecs = mdpspecs['defaults'] if 'defaults' in mdpspecs else ()
			for mdpname in [i for i in mdpspecs if re.match('.+\.mdp$',i)]:
				print 'generating mdp: '+mdpname
				#---for each mdp file the specs defined in the dictionary help select between redundant options
				specs = mdpspecs[mdpname]+defspecs
				routes = list(mapdict(delve(mdpdict,*topkeys)))
				#---any incomplete routes (partials) imply inclusion of all routes that match
				partials = [spec for spec in specs if spec not in routes]
				valids = [r for r in routes if 
					(any([r[:len(p)]==p for p in partials]) and #---under a partial
					r[:-1] not in [s[:-1] for s in specs]) or #---no conflicts
					r in specs] #---explicit route in specs
				keycollect = [delve(mdpdict,*(topkeys+r)).keys() for r in valids]
				if redundant([i for j in keycollect for i in j]): 
					raise Exception('repeated mdp keys: '+str(valids))
				with open(self.rootdir+mdpname,'w') as fp:
					for v in valids:
						for key,val in delve(mdpdict,*(topkeys+v)).items():
							fp.write(str(key)+' = '+str(val)+'\n')
		
		if 0:
			#--topkeys is the root node for our parameters in mdpdict
			topkeys = mdpspecs['top'] if 'top' in mdpspecs else ()
			#---defspecs set any default paths on the mdpdict tree which we use for all mdp files
			defspecs = mdpspecs['defaults'] if 'defaults' in mdpspecs else ()
			for mdpname in [i for i in mdpspecs if re.match('.+\.mdp$',i)]:
				print 'generating mdp: '+mdpname
				#---for each mdp file the specs defined in the dictionary help select between redundant options
				specs = mdpspecs[mdpname]+defspecs
				routes = list(mapdict(delve(mdpdict,*topkeys)))
				valids = [r for r in routes if r in specs or r[:-1] not in [s[:-1] for s in specs]]
				keycollect = [delve(mdpdict,*(topkeys+r)).keys() for r in valids]
				if redundant([i for j in keycollect for i in j]): raise Exception('repeated mdp keys: '+str(valids))
				with open(self.rootdir+mdpname,'w') as fp:
					for v in valids:
						for key,val in delve(mdpdict,*(topkeys+v)).items():
							fp.write(str(key)+' = '+str(val)+'\n')
	
		if 0:
			mdpspecs = {}
			execfile(os.path.abspath(os.path.expanduser('./inputs/input-specs-mdp.dat')),mdpspecs)
			mdpdict = mdpspecs['mdpdict']
			for mdpname,specs in self.params['bilayer_construction_settings'][self.simscale]['mdp'].items():
				print 'generating MDP: '+mdpname
				print specs
				for key in specs:
					print key
					print specs[key]
			bp()

class Multiply(AMXSimulation):

	"""
	A generic class which continues a larger, periodic replicate of a simulation.
	"""
	
	def __init__(self,rootdir=None,previous_dir=None,nx=1,ny=1,nz=1,**kwargs):

		#---set paths for the previous and current steps
		self.rootdir = os.path.abspath(os.path.expanduser(rootdir))+'/'
		self.prevdir = os.path.abspath(os.path.expanduser(previous_dir))+'/'
		sys.stdout = tee(open(self.rootdir+'log-script-master','a',1))
		sys.stderr = tee(open(self.rootdir+'log-script-master','a',1),error=True)
		print 'MULTIPLYING THE SIMULATION BOX'

		#---search backwards to find the most recent topology
		startstep,oldsteps = chain_steps()
		for ostep in oldsteps[::-1]:
			for root,dirnames,filenames in os.walk(ostep): break
			if 'system.top' in filenames: break
		for fn in ['system.top','martini.ff','charmm36.ff','lipids-tops','input-md-in.mdp']:
			if os.path.isfile(self.rootdir+'../'+ostep+'/'+fn) or \
				os.path.isdir(self.rootdir+'../'+ostep+'/'+fn):
				copy(self.rootdir+'../'+ostep+'/'+fn,self.rootdir+fn)
		copy(self.rootdir+'system.top',self.rootdir+'system-small.top')
		with open(self.rootdir+'system-small.top','r') as fp: small = fp.readlines()
		startmol = [ii for ii,i in enumerate(small) if re.match('^(\s+)?\[\s?molecules\s?\]',i)][0]
		
		#---multiply the topology
		nx = int(nx) if nx != '' else 1
		ny = int(ny) if ny != '' else 1
		nz = int(nz) if nz != '' else 1
		newtop = []
		for s in small[startmol+1:]:
			if re.match('^[A-Z,a-z,_,+,\-,0-9]+\s+[0-9]+$',s.strip('\n')):
				newtop.append([s.split()[0],str(int(s.split()[1])*(nx*ny))])
		with open(self.rootdir+'system.top','w') as fp:
			for i in small[:startmol+1]: fp.write(i)
			for n in newtop: fp.write(n[0]+' '+n[1]+'\n')
		newtop = [[str(i[0]),int(i[1])] for i in newtop]
		self.ion_residue_names = ['W','SOL','MG','NA','CL','Cal']
		self.lnames = [i[0] for i in newtop]
				
		#---get the last frame using the second-to-last directory now that we made a multiply directory
		self.smallstruct = latestcheck(oldsteps[-2])[0][:-4]+'.gro'
		lastframe(rootdir=oldsteps[-2],prefix=self.smallstruct[:-4],gmxpaths=gmxpaths)
		copy(oldsteps[-2]+'/'+self.smallstruct,self.rootdir+self.smallstruct)
			
		#---multiply the structure
		cmd = [gmxpaths['genconf'],
			'-f '+self.smallstruct,
			'-o system-input-unsorted.gro',
			'-nbox '+' '.join([str(i) for i in [nx,ny,nz]])]
		call(cmd,logfile='log-genconf',cwd=self.rootdir)
		
		#---reorder the gro file
		with open(self.rootdir+'system-input-unsorted.gro','r') as fp: lines = fp.readlines()
		reord = lines[:2]
		for res in newtop:
			for line in lines:
				if line[5:10].strip() == res[0] or line[10:15].strip() == res[0]: reord.append(line)
		reord.append(lines[-1])
		
		with open(self.rootdir+'system-input-unordered.gro','w') as fp:
			for line in reord: fp.write(line)

		#---multiply the structure
		cmd = [gmxpaths['editconf'],
			'-f system-input-unordered.gro',
			'-o system-input.gro',
			'-resnr 1']
		call(cmd,logfile='log-editconf-resnr',cwd=self.rootdir)
		self.grouping(grouptype='bilayer',startstruct='system-input.gro')
		os.remove(self.rootdir+self.smallstruct)

