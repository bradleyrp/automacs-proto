#!/usr/bin/python -B

#---starting structure
input_filename = 'sources/aamd-protein-structs/3ala.pdb'

#---controller sets the scale automatically
simscale = 'cgmd'

#---settings for protein.py
protein_construction_settings = {
	'aamd':{
		'positive_counterion_name':'NA',
		'cube':False,
		'ion_strength':0.15,
		'solvent_structure':'spc216.gro',
		'align_x':True,
		'water_model':'tip3p',
		'negative_counterion_name':'CL',
		'histype':'d',
		'sol_name':'SOL',
		'wbuffer':1.5,
		'force_field':'charmm27',
		'force_field_local':None,
		'top_includes':[
			'charmm27.ff/forcefield.itp',
			'protein.itp',
			'charmm27.ff/tip3p.itp',
			'charmm27.ff/ions.itp',
			],
		'mdp':{
			'group':'aamd',
			'input-em-steep-in.mdp':['minimize'],
			'input-em-cg-in.mdp':['minimize',{'integrator':'cg'}],
			},
		'mdp-equil':{
			'group':'aamd',
			'input-md-in.mdp':None,
			'input-md-nvt-eq-in.mdp':['nvt-protein'],
			'input-md-npt-eq-in.mdp':['npt-protein'],
			},
		},
	'cgmd':{
		'positive_counterion_name':'NA+',
		'cube':False,
		'protein_water_gap':3,
		'ion_strength':0.15,
		'solvent_structure':'solvate-water.gro',
		'negative_counterion_name':'CL-',
		'sol_name':'W',
		'wbuffer':3,
		'top_includes':[
			'martini.ff/martini-v2.2.itp',
			'martini.ff/martini-v2.0-lipids.itp',
			'martini.ff/martini-v2.2-aminoacids.itp',
			'martini.ff/martini-v2.0-ions.itp'
			]
		}
	}