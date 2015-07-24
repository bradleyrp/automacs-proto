#!/usr/bin/python -B

#---COMPOSITIONS
nlipids0 = 500
complist0 = [4,1,0]
lipidnames0 = ['DOPC','DOPS','CHOL']

#---ASYMMETRIC COMPOSITION (omit for symmetric)
nlipids1 = 500
complist1 = [4,1,0]
lipidnames1 = ['DOPC','DOPS','CHOL']

#---SIMULATION RESOLUTION
simscale = 'aamd'

#---CONSTRUCTION SETTINGS
bilayer_construction_settings = {
	'aamd': {
		'lbuffer':0.1,
		'tbuffer':0.1,
		'zbuffer':1.5,
		'square':True,
		'solvent_thickness':10,
		'lipid_water_buffer':0.1,
		'ion_strength':0.150,
		'sol_name':'SOL',
		'negative_ion_name':'CL','nq':-1,
		'positive_ion_name':'NA','pq':1,
		'lipid_spacing':0.2,
		'vacuum_packing_sequence':
			['input-md-vacuum-p1-in.mdp',
			'input-md-vacuum-p2-in.mdp'],
		'water_conf':'spc216.gro',
		'ion_residue_names':['MG','NA','CL','Cal'],
		'concentration_calc':['simple','lipid_exempt'][0],
		'top_includes':[],
		},
	'cgmd': {
		'lbuffer':0.1,
		'tbuffer':0.1,
		'zbuffer':4.0,
		'square':True,
		'solvent_thickness':10,
		'lipid_water_buffer':0.1,
		'ion_strength':0.150,
		'sol_name':'W',
		'negative_ion_name':'CL-','nq':-1,
		'positive_ion_name':'NA+','pq':1,
		'salt':['NA+','CL-'],
		'lipid_spacing':0.2,
		'vacuum_packing_sequence':
			['input-md-vacuum.mdp',
			'input-md-vacuum2.mdp',
			'input-md-vacuum3.mdp'],
		'ion_residue_names':['NA+','CL-','ION'],
		'concentration_calc':['simple','lipid_exempt'][0],
		'water_conf':'structures/martini-water.gro',
		'top_includes':[
			'martini.ff/martini-v2.2.itp',
			'martini.ff/martini-v2.0-lipids.itp',
			'lipids-tops/PIP2.itp',
			'martini.ff/martini-v2.2-aminoacids.itp',
			'martini.ff/martini-v2.0-ions.itp',
			'martini.ff/martini-v2.0-cholesterol.itp',
			],
		'mdp':{
			'group':'cgmd',
			'input-em-steep-in.mdp':['minimize'],
			'input-em-cg-in.mdp':['minimize',{'integrator':'cg'}],
			'input-em-steep-posre-in.mdp':['minimize',{'restrain':'posre'}],
			'input-em-cg-posre-in.mdp':['minimize',{'integrator':'cg'},{'restrain':'posre'}],
			'input-md-vacuum.mdp':['vacuum-packing',{'ref_p':'100.0 1.0','dt':0.001}],
			'input-md-vacuum2.mdp':['vacuum-packing',{'ref_p':'10.0 1.0','dt':0.002}],
			'input-md-vacuum3.mdp':['vacuum-packing'],
			},
		'mdp-equil':{
			'group':'cgmd',
			'input-md-in.mdp':None,
			'input-md-npt-eq-in.mdp':['npt-bilayer'],
			},
		},
	}
