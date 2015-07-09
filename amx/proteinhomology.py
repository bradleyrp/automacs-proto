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
import glob

#---locate gmxpaths
if os.path.isfile('gmxpaths.conf'): gmxpaths_path = './'
else: gmxpaths_path = '../'

#---load the gromacs paths
gmxpaths = dict([[i.split()[0],' '.join(i.split()[1:])] 
	for i in open(gmxpaths_path+'gmxpaths.conf','r') 
	if i.strip()[0] != '#'])
	
dna_mapping = {
	'UUU':'F','UUC':'F','UUA':'L','UUG':'L','UCU':'S','UCC':'s','UCA':'S','UCG':'S','UAU':'Y','UAC':'Y',
	'UAA':'STOP','UAG':'STOP','UGU':'C','UGC':'C','UGA':'STOP','UGG':'W','CUU':'L','CUC':'L','CUA':'L',
	'CUG':'L','CCU':'P','CCC':'P','CCA':'P','CCG':'P','CAU':'H','CAC':'H','CAA':'Q','CAG':'Q','CGU':'R',
	'CGC':'R','CGA':'R','CGG':'R','AUU':'I','AUC':'I','AUA':'I','AUG':'M','ACU':'T','ACC':'T','ACA':'T',
	'ACG':'T','AAU':'N','AAC':'N','AAA':'K','AAG':'K','AGU':'S','AGC':'S','AGA':'R','AGG':'R','GUU':'V',
	'GUC':'V','GUA':'V','GUG':'V','GCU':'A','GCC':'A','GCA':'A','GCG':'A','GAU':'D','GAC':'D','GAA':'E',
	'GAG':'E','GGU':'G','GGC':'G','GGA':'G','GGG':'G',} 
	
aacodemap = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
	'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER':'*',
	'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA':'X'}

class ProteinHomology:
	def __init__(self,rootdir=None,**kwargs):
	
		#---set paths for the previous and current steps
		self.rootdir = os.path.abspath(os.path.expanduser(rootdir))+'/'
		self.rootdirrel = rootdir+'/'
		
		#---set sources
		self.sources_dir = os.path.abspath(os.path.expanduser('./sources/'))+'/'
		
		#---use default inputs file
		self.inputs_file = os.path.abspath(os.path.expanduser('./inputs/input_specs_homology.py'))
		
		#---use default path to store structure files
		self.struct_dir = os.path.abspath(os.path.expanduser(
			self.sources_dir+'scripts/aamd-protein-structs'))+'/'
		
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
		copy(self.sources_dir+'scripts/aamd-protein-homology/*',self.rootdir)

		#---run the homology building procedure
		self.build()
		
	def build(self):
	
		print 'BUILDING HOMOLOGY MODELS'
		if self.procedure != 'mutator': self.get_targets()
		self.get_templates()
		if self.procedure == 'single': self.build_model_single()
		elif self.procedure == 'multi': self.build_model_multi()
		elif self.procedure == 'mutator': self.build_model_mutator()
		
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
			elif key == 'textfile_rna':
				with open(target_ins[key],'r') as fp: targs = fp.readlines()
				for t in targs:
					if re.match('^[a-z,A-Z,0-9,_].+\s*:\s*[A-Z,a-z].+$',t):
						self.target.append(list([i.strip() for i in t.split(':')]))
						rnaseq = self.target[-1][1]
						#---extra substitutions for later
						if 'regex_subs' in self.settings.keys():
							for regex in self.settings['regex_subs']:
								rnaseq = re.sub(regex[0],regex[1],rnaseq)
						rnaseq = rnaseq.upper()
						rnaseq = re.sub('T','U',rnaseq)
						aminoseq = ''.join([dna_mapping[i] for i in [rnaseq[i:i+3] 
							for i in range(0,len(rnaseq),3)]])
						self.target[-1][1] = re.sub('T','U',aminoseq)
						self.target[-1] = tuple(self.target[-1])
			else: raise Exception('except: unclear target type')

	def get_templates(self):

		if not os.path.isdir('./repo'): os.mkdir('./repo')
		temps = self.settings['template']
		#---ensure that the template object is always in a list
		if len(temps) == 2 and type(temps[0])==str and type(temps[1])==str: temps = [temps]
		self.template = []
		for t in temps:
			print 'retrieving '+str(t[0])
			#---check if in repo and move
			if not os.path.isfile(self.rootdir+t[0]+'.pdb') and os.path.isfile('./repo/'+t[0]+'.pdb'):
				copy('./repo/'+t[0]+'.pdb',self.rootdir+t[0]+'.pdb')
				#---fasta retrieval is deprecated
				if 0: copy('./repo/'+t[0]+'.fasta',self.rootdir+t[0]+'.fasta')
			elif not os.path.isfile(self.rootdir+t[0]+'.pdb'):
				response = urllib2.urlopen('http://www.rcsb.org/pdb/files/'+t[0]+'.pdb')
				pdbfile = response.read()
				with open(self.rootdir+t[0]+'.pdb','w') as fp: fp.write(pdbfile)
				copy(self.rootdir+t[0]+'.pdb','./repo/'+t[0]+'.pdb')
				#---always get the FASTA sequence whenever you get a structure
				#---note currently deprecated
				if 0:
					response = urllib2.urlopen(
						'http://www.rcsb.org/pdb/files/fasta.txt?structureIdList='+t[0])
					fastafile = response.read()
					with open(self.rootdir+t[0]+'.fasta','w') as fp: fp.write(fastafile)
					copy(self.rootdir+t[0]+'.fasta','./repo/'+t[0]+'.fasta')
			self.template.append(t)
			
	def build_model_single(self,targi=0,tempi=0,batchdir_override=None):
	
		if len(self.template) != 1: raise Exception('except: needs only one template '+str(self.template))
		if 'target_name' in self.settings.keys():
			if type(self.settings['target_name']) == str: 
				self.settings['target_name'] = [self.settings['target_name']]
		
		if len(self.target) == 1: targi_inds = [0]
		elif len(self.target)>1 and 'target_name' in self.settings.keys() and \
			len(self.settings['target_name']) == 1:
			targi_inds = [[i[0] for i in self.target].index(self.settings['target_name'][0])]
		elif len(self.target)>1 and 'target_name' in self.settings.keys() and \
			len(self.settings['target_name']) > 1:
			targi_inds = [[i[0] for i in self.target].index(j) for j in self.settings['target_name']]
		elif len(self.target)>1 and 'target_name' not in self.settings.keys():
			print 'multiple targets selected so this will be a batch operation'
			targi_inds = range(len(self.target))
		else: raise Exception('except: target/target_name mismatch')
		
		for targi in range(len(targi_inds)):
			if len(targi_inds)>1:
				print 'BATCH MODEL GENERATION, MODEL No. = '+str(targi)
				batchdir = self.rootdir+'model-v'+('%04d'%(targi+1))+'-'+self.target[targi][0]+'/'
				batchdirrel = self.rootdirrel+'model-v'+('%04d'%(targi+1))+'-'+self.target[targi][0]+'/'
				os.mkdir(batchdir)
				copy(self.sources_dir+'scripts/aamd-protein-homology/*',batchdir)
				copy(self.rootdir+'*.pdb',batchdir)
			elif batchdir_override != None:
				batchdir = self.rootdir+batchdir_override+'/'
				batchdirrel = self.rootdirrel+batchdir_override+'/'
				os.mkdir(batchdir)
				copy(self.sources_dir+'scripts/aamd-protein-homology/*',batchdir)
				copy(self.rootdir+'*.pdb',batchdir)
			else: batchdir,batchdirrel = self.rootdir,self.rootdirrel
		
			print 'preparing modeller scripts'
			#---variables passed to modeller via settings-homology.py
			vars_to_modeller = {
				'template_struct':self.template[0][0],
				'template_struct_chain':self.template[0][1],
				'target_seq':self.target[targi][0],
				'n_models':self.settings['n_models'],
				}
	
			#---write a settings file for the modeller script
			with open(batchdir+'settings-homology.py','w') as fp:
				fp.write('#!/usr/bin/python\n\n')
				for var in vars_to_modeller.keys():
					val = '\''+str(vars_to_modeller[var])+'\'' \
						if type(vars_to_modeller[var]) == str else vars_to_modeller[var]
					fp.write(var+' = '+str(val)+'\n')
			
			#---write an ali file with the target
			fasta_linelen = 50
			with open(batchdir+self.target[targi][0]+'.ali','w') as fp:
				fp.write('>P1;'+self.target[targi][0]+'\n')
				fp.write('sequence:'+self.target[targi][0]+':::::::0.00:0.00\n')
				seq = self.target[targi][1]
				chopped = [seq[j*fasta_linelen:(j+1)*fasta_linelen] for j in range(len(seq)/fasta_linelen+1)]
				chopped = [i for i in chopped if len(i) > 0]
				for i,seg in enumerate(chopped): fp.write(seg+('\n' if i < len(chopped)-1 else '*\n'))
		
			print 'running modeller'
			cmd = [gmxpaths['modeller'],'script-single.py']
			call(cmd,logfile='log-modeller-script-single',cwd=batchdir)
		
			#---make a view script
			with open(batchdir+'script-vmd-view.tcl','w') as fp:
				fp.write('#!/usr/bin/env tclsh\n\n#---execute via "vmd -e script-vmd-view.tcl"\n\n')			
				fp.write('mol default style {NewCartoon 0.300000 6.000000 4.100000 0}\n')
				fp.write('mol default color {Structure}\n')
				for fn in glob.glob(batchdirrel+self.target[targi][0]+'.*.pdb'):
					fp.write('mol new '+os.path.basename(fn)+'\n')
				fp.write('for {set i 0} { $i <= '+str(self.settings['n_models'])+'} { incr i 1} {\n')
				fp.write('\tset reference_sel  [atomselect 1 "backbone"]\n')
				fp.write('\tset comparison_sel [atomselect $i "backbone"]\n')
				fp.write('\tset transformation_mat [measure fit $comparison_sel $reference_sel]\n')
				fp.write('\tset move_sel [atomselect $i "all"]\n')
				fp.write('\t$move_sel move $transformation_mat\n}\n')
		
		#---write comparison script here
				
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


	def build_model_mutator(self):

		print "beginning mutation procedure"
		#---deprecated use of fasta file for getting the sequence
		if 0:
			with open(self.rootdir+self.template[0][0]+'.fasta','r') as fp: lines = fp.readlines()
			seqstarts = [li for li,l in enumerate(lines) if re.match('^>',l)]+[len(lines)]
			pdbcodes = [re.findall('^>([A-Z,0-9]{4})\:([A-Z]{1})\|',lines[i])[0] for i in seqstarts[:-1]]
			seqs = [''.join([lines[l].strip('\n') for l in range(seqstarts[i-1]+1,seqstarts[i])]) 
				for i in range(1,len(seqstarts))]
			template = self.template[0]
			sequence = seqs[pdbcodes.index(template)]
			#---note that this is a hack but we get the first residue index from the PDB file directly
			with open(self.rootdir+template[0]+'.pdb','r') as fp: lines = fp.readlines()
			startres = int([l for l in lines if re.match('^ATOM',l)][0].split()[5])
		with open(self.rootdir+self.template[0][0]+'.pdb','r') as fp: lines = fp.readlines()
		seqresli = [li for li,l in enumerate(lines) if re.match('^SEQRES\s+',l)]
		seqraw = [re.findall('^SEQRES\s+[0-9]+\s+([A-Z])\s+[0-9]+\s+(.+)',lines[li])[0] for li in seqresli]
		sequence = ''.join([''.join([aacodemap[j] for j in i[1].split()]) 
			for i in seqraw if i[0] == self.template[0][1]])
		#---handle missing residues
		missingli = [re.findall('^REMARK\s+([0-9]+)\sMISSING RESIDUES',l)[0] for li,l in enumerate(lines) 
			if re.match('^REMARK\s+([0-9]+)\sMISSING RESIDUES',l)][0]
		startres = int([re.findall('^REMARK\s+'+missingli+'\s+[A-Z]{3}\s+[A-Z]\s+([0-9]+)',l)[0]
			for li,l in enumerate(lines) 
			if re.match('^REMARK\s+'+missingli+'\s+[A-Z]{3}\s+[A-Z]\s+[0-9]+',l)][0])
		self.target = []
		for mi,mut in enumerate(self.settings['mutations']):
			sequence_mut = list(sequence)
			if sequence[mut[1]-startres] != mut[0]: raise Exception('wrong starting residue')
			else: sequence_mut[mut[1]-startres] = mut[2]
			sequence_mut = ''.join(sequence_mut)
			print 'template sequence = '+sequence
			print 'mutated sequence = '+sequence_mut
			self.target.append(['mutation'+str(mi),sequence_mut])
		for mi,mut in enumerate(self.settings['mutations']):
			print 'building homology model for mutation '+str(mi)
			self.settings['target_name'] = self.target[mi][0]
			batchdir = 'model-v'+('%04d'%(mi+1))+'-'+''.join(self.template[0])+\
				'_'+''.join([str(j) for j in mut])+'/'
			self.build_model_single(targi=mi,tempi=0,batchdir_override=batchdir)

		
