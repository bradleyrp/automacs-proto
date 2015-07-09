from modeller import *
from modeller.automodel import *
import sys

#---SETTINGS
#-------------------------------------------------------------------------------------------------------------

execfile('settings-homology.py')
doalign2d = True

#---2D ALIGNMENT
#-------------------------------------------------------------------------------------------------------------

if doalign2d:
	env = environ()
	aln = alignment(env)
	mdl = model(env, \
		file=template_struct, \
		model_segment=('FIRST:'+template_struct_chain, 'LAST:'+template_struct_chain))
	aln.append_model(mdl, \
		align_codes=template_struct+template_struct_chain, \
		atom_files=template_struct+'.pdb')
	aln.append(file=target_seq+'.ali', align_codes=target_seq)
	aln.align2d()
	aln.write(file='align2d.ali', alignment_format='PIR')
	aln.write(file='align2d.pap', alignment_format='PAP')
	afile = 'align2d.ali'
else:
	env = environ()
	aln = alignment(env)
	afile = 'align2d-custom.ali'

#---BUILD MODELS
#-------------------------------------------------------------------------------------------------------------

a = automodel(env, \
	alnfile=afile, \
	knowns=template_struct+template_struct_chain, \
	sequence=target_seq, \
	assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1
a.ending_model = n_models
a.make()

