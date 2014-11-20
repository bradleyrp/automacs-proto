#!/usr/bin/python

from tools import call,checkout,tee,copy,chain_steps,latestcheck,lastframe
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

'''
Parent class for any GROMACS simulation. Comprised of common functions necessary for building most systems.
'''

class AMXSimulation:

	def get_box_vectors(self,logfile,new=True):
		'''Read box vectors from editconf output stored in a log file.'''
		if new:
			boxdims = [float(checkout(["awk","'/new box vectors/ {print $"+\
				str(i+5)+"}'",logfile],cwd=self.rootdir).strip()) for i in range(3)]
		else:
			boxdims = [float(checkout(["awk","'/    box vectors/ {print $"+\
				str(i+4)+"}'",logfile],cwd=self.rootdir).strip()) for i in range(3)]
		return boxdims
		
	def write_topology_bilayer(self,topname,water_only=False):
		'''
		Write the bilayer topology file from residue lists stored in ``self.lnames`` and ``self.complist``.
		'''
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

	def write_topology_protein(self,topname):
		'''
		Write a protein+water topology file.
		'''
		print "printing protein topology file"
		if self.simscale == 'aamd':
			fp = open(self.rootdir+topname,'w')
			fp.write('#include "charmm27.ff/forcefield.itp"\n')
			fp.write('#include "protein.itp"\n')
			fp.write('; Include water topology\n')
			fp.write('#include "charmm27.ff/tip3p.itp"\n')
			fp.write('#ifdef POSRES_WATER\n')
			fp.write('; Position restraint for each water oxygen\n')
			fp.write('[ position_restraints ]\n')
			fp.write(';  i funct       fcx        fcy        fcz\n')
			fp.write('   1    1       1000       1000       1000\n')
			fp.write('#endif\n')
			fp.write('; Include topology for ions\n')
			fp.write('#include "charmm27.ff/ions.itp"\n\n')
			fp.write('[ system ]\nPROTEIN+WATER\n\n')
			fp.write('[ molecules ]\n')
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
			fp.write('#include "martini.ff/martini-v2.2.itp"'+'\n')
			fp.write('#include "martini.ff/martini-v2.0-lipids.itp"'+'\n')
			for lipid_itp in self.itp_lipid: fp.write('#include "lipids-tops/'+lipid_itp+'.itp"'+'\n')
			fp.write('#include "martini.ff/martini-v2.2-aminoacids.itp"'+'\n')
			fp.write('#include "martini.ff/martini-v2.0-ions.itp"'+'\n')
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
		'''
		Write the ``system-groups.ndx`` file for further simulation, which requires that the water and 
		lipids and possibly proteins are coupled to their own temperature baths.
		'''
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
				' | r '.join([l for l in self.lnames if l not in self.ion_residue_names]+['ION'])+\
				'\n'+'r '+' | r '.join([l for l in self.lnames if l in self.ion_residue_names])+\
				'\nname 1 LIPIDS\nname 2 SOLV\ndel 0\nq\n'
			#---write a groups file suitable for bilayers with separate lipid and solvent groups
			cmd = [gmxpaths['make_ndx'],
				'-f '+startstruct,
				'-o system-groups.ndx']
			call(cmd,logfile='log-make-ndx-groups',cwd=self.rootdir,inpipe=inpipe)
		else: raise Exception('except: incomprehensible group type')

class Multiply(AMXSimulation):

	'''A generic class which continues a larger, periodic replicate of a simulation.'''
	
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
			if re.match('^[A-Z,a-z,_,+,-]+\s+[0-9]+$',s.strip('\n')):
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
		
		
		

