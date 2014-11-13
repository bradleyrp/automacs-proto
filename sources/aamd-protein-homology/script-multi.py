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
	log.verbose()
	env = environ()
	env.io.atom_files_directory = './:../atom_files/'
	aln = alignment(env)
	for code in pdblist:
		pdb = code[0]
		chain = code[1]
		mdl = model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
		aln.append_model(mdl, atom_files=pdb, align_codes=pdb+chain)
	for (weights, write_fit, whole) in \
		(((1., 0., 0., 0., 1., 0.), False, True),
		((1., 0.5, 1., 1., 1., 0.), False, True),
		((1., 1., 1., 1., 1., 0.), True, False)):
			aln.salign(rms_cutoff=3.5, 
			normalize_pp_scores=False,
			rr_file='$(LIB)/as1.sim.mat', 
			overhang=30,
			gap_penalties_1d=(-450, -50),
			gap_penalties_3d=(0, 3), 
			gap_gap_score=0, 
			gap_residue_score=0,
			dendrogram_file='dendrogram.tree',
			alignment_type='tree',feature_weights=weights,
			improve_alignment=True, fit=True, write_fit=write_fit,
			write_whole_pdb=whole, 
			output='ALIGNMENT QUALITY')	
	aln.write(file='alignment.pap', alignment_format='PAP')
	aln.write(file='alignment.ali', alignment_format='PIR')
	aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
		rr_file='$(LIB)/as1.sim.mat', overhang=30,
		gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
		gap_gap_score=0, gap_residue_score=0,
		alignment_type='progressive', feature_weights=[0]*6,
		improve_alignment=False, fit=False, write_fit=True,
		write_whole_pdb=False, output='QUALITY')
	log.verbose()
	env2 = environ()
	env2.libs.topology.read(file='$(LIB)/top_heav.lib')
	aln2 = alignment(env2)
	aln2.append(file='alignment.ali', align_codes='all')
	aln_block = len(aln2)
	aln2.append(file=target_seq+'.ali', align_codes=target_seq)
	aln2.salign(output='', max_gap_length=20,
		gap_function=True,   
		alignment_type='PAIRWISE', align_block=aln_block,
		feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
		gap_penalties_1d=(-450, 0),
		gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
		similarity_flag=True)
	aln2.write(file='align-mult.ali', alignment_format='PIR')
	aln2.write(file='align-mult.pap', alignment_format='PAP')
	afile = 'align-mult.ali'
else:
	afile = 'align-mult-custom.ali'

#---BUILD MODELS
#-------------------------------------------------------------------------------------------------------------

pdblist2 = []
for item in pdblist:
	pdblist2.append(item[0]+item[1])
env3 = environ()
a = automodel(env3, alnfile=afile,
	knowns=pdblist2, sequence=target_seq)
a.starting_model = 1
a.ending_model = n_models
a.make()	

#---------------------------------------------------------------------------------------------------
