#!/usr/bin/python

#---AUTOMACS (AMX) "PROTEIN HOMOLOGY" MODULE
#---created 2014.11.09 by ryan bradley

#---development features
import os,subprocess
if os.path.isfile('/etc/pythonstart'):
	execfile('/etc/pythonstart')

#---imports
import os,sys,time,datetime,re
import urllib2
from tools import call,checkout,tee,copy
import amxsim

#---locate gmxpaths
if os.path.isfile('gmxpaths.conf'): gmxpaths_path = './'
else: gmxpaths_path = '../'

#---load the gromacs paths
gmxpaths = dict([[i.split()[0],' '.join(i.split()[1:])] 
	for i in open(gmxpaths_path+'gmxpaths.conf','r') 
	if i.strip()[0] != '#'])

class ProteinHomology:
	def __init__(self,rootdir=None,**kwargs):
	
		#---set paths for the previous and current steps
		self.rootdir = os.path.abspath(os.path.expanduser(rootdir))+'/'
		
		#---set sources
		self.sources_dir = os.path.abspath(os.path.expanduser('./sources/'))+'/'
		
		#---use default inputs file
		self.inputs_file = os.path.abspath(os.path.expanduser('./inputs/input-specs-homology.dat'))
		
		#---use default path to store structure files
		self.struct_dir = os.path.abspath(os.path.expanduser(self.sources_dir+'aamd-protein-structs'))+'/'
		
		#---protein construction settings are pulled from the inputs_file
		self.params = dict()
		execfile(self.inputs_file,self.params)
		self.procedure = self.params['procedure']
		self.settings = self.params['homology_construction_settings'][self.procedure]
		
		#---make root directory
		if not os.path.isdir(self.rootdir): os.mkdir(rootdir)

		#---start writing to log-script-master		
		sys.stdout = tee(open(self.rootdir+'log-script-master','a',1))
		sys.stderr = tee(open(self.rootdir+'log-script-master','a',1),error=True)
		
		#---copy modeller scripts
		copy(self.sources_dir+'aamd-protein-homology/*',self.rootdir)

		#---run the homology building procedure
		self.build()
		
	def build(self):
	
		print 'BUILDING HOMOLOGY MODELS'
		self.get_targets()
		self.get_templates()
		if self.procedure == 'single': self.build_model_single()
		elif self.procedure == 'multi': self.build_model_multi()
		
	def get_targets(self):
	
		self.target = []
		target_ins = self.settings['target']
		for key in target_ins.keys():
			if key == 'raw':
				self.target.append(target_ins[key])
			elif key == 'textfile':
				with open(target_ins[key],'r') as fp: targs = fp.readlines()
				for t in targs:
					if re.match('^[a-z,A-Z,_].+\s*:\s*[A-Z].+$',t):
						self.target.append(tuple([i.strip() for i in t.split(':')]))
			else: raise Exception('except: unclear target type')
		
	def get_templates(self):

		if not os.path.isdir('./repo'): os.mkdir('./repo')
		temps = self.settings['template']
		self.template = []
		for t in temps:
			print 'retrieving '+str(t[0])
			#---check if in repo and move
			if not os.path.isfile(self.rootdir+t[0]+'.pdb') and os.path.isfile('./repo/'+t[0]+'.pdb'):
				copy('./repo/'+t[0]+'.pdb',self.rootdir+t[0]+'.pdb')			
			elif not os.path.isfile(self.rootdir+t[0]+'.pdb'):
				response = urllib2.urlopen('http://www.rcsb.org/pdb/files/'+t[0]+'.pdb')
				pdbfile = response.read()
				with open(self.rootdir+t[0]+'.pdb','w') as fp: 
					fp.write(pdbfile)
				copy(self.rootdir+t[0]+'.pdb','./repo/'+t[0]+'.pdb')
			self.template.append(t)
			
	def build_model_single(self):
	
		if len(self.template) != 1: raise Exception('except: needs only one template '+str(self.template))
		if len(self.target) != 1: raise Exception('except: needs only one target '+str(self.target))
	
		print 'preparing modeller scripts'
		#---variables passed to modeller via settings-homology.py
		vars_to_modeller = {
			'template_struct':self.template[0][0],
			'template_struct_chain':self.template[0][1],
			'target_seq':self.target[0][0],
			'n_models':self.settings['n_models'],
			}
	
		#---write a settings file for the modeller script
		with open(self.rootdir+'settings-homology.py','w') as fp:
			fp.write('#!/usr/bin/python\n\n')
			for var in vars_to_modeller.keys():
				val = '\''+str(vars_to_modeller[var])+'\'' \
					if type(vars_to_modeller[var]) == str else vars_to_modeller[var]
				fp.write(var+' = '+str(val)+'\n')
			
		#---write an ali file with the target
		fasta_linelen = 50
		with open(self.rootdir+self.target[0][0]+'.ali','w') as fp:
			fp.write('>P1;'+self.target[0][0]+'\n')
			fp.write('sequence:'+self.target[0][0]+':::::::0.00:0.00\n')
			seq = self.target[0][1]
			chopped = [seq[j*fasta_linelen:(j+1)*fasta_linelen] for j in range(len(seq)/fasta_linelen+1)]
			chopped = [i for i in chopped if len(i) > 0]
			for i,seg in enumerate(chopped): fp.write(seg+('\n' if i < len(chopped)-1 else '*\n'))
		
		print 'running modeller'
		cmd = [gmxpaths['modeller'],'script-single.py']
		call(cmd,logfile='log-modeller-script-single',cwd=self.rootdir)
		
	def build_model_multi(self):
	
		if len(self.template) < 1: raise Exception('except: needs multiple templates '+str(self.template))
		if len(self.target) != 1: raise Exception('except: needs only one target '+str(self.template))
	
		print 'preparing modeller scripts'
		#---variables passed to modeller via settings-homology.py
		vars_to_modeller = {
			'pdblist':self.template,
			'target_seq':self.target[0][0],
			'n_models':self.settings['n_models'],
			}
	
		#---write a settings file for the modeller script
		with open(self.rootdir+'settings-homology.py','w') as fp:
			fp.write('#!/usr/bin/python\n\n')
			for var in vars_to_modeller.keys():
				val = '\''+str(vars_to_modeller[var])+'\'' \
					if type(vars_to_modeller[var]) == str else vars_to_modeller[var]
				fp.write(var+' = '+str(val)+'\n')
			
		#---write an ali file with the target
		fasta_linelen = 50
		with open(self.rootdir+self.target[0][0]+'.ali','w') as fp:
			fp.write('>P1;'+self.target[0][0]+'\n')
			fp.write('sequence:'+self.target[0][0]+':::::::0.00:0.00\n')
			seq = self.target[0][1]
			chopped = [seq[j*fasta_linelen:(j+1)*fasta_linelen] for j in range(len(seq)/fasta_linelen+1)]
			chopped = [i for i in chopped if len(i) > 0]
			for i,seg in enumerate(chopped): fp.write(seg+('\n' if i < len(chopped)-1 else '*\n'))
		
		print 'running modeller'
		cmd = [gmxpaths['modeller'],'script-multi.py']
		call(cmd,logfile='log-modeller-script-multi',cwd=self.rootdir)

