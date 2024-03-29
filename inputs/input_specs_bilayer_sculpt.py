#!/usr/bin/python

#---SIMULATION RESOLUTION
simscale = 'aamd'

#---CONSTRUCTION SETTINGS

bilayer_construction_settings = {
	'cgmd': {
		'sol_name':'W',
		'simuluxe_script':'~/worker/simuluxe/project-bilayer/construct-shapes/script-sculptor.py',
		'simuluxe_script_restrain':'~/worker/simuluxe/project-bilayer/construct-shapes/script-fixer.py',
		'negative_ion_name':'CL-','nq':-1,
		'positive_ion_name':'NA+','pq':1,
		'solvent_structure':'solvate-water.gro',
		'ion_strength':0.150,
		'ion_residue_names':['NA+','CL-','ION'],
		'concentration_calc':['simple','lipid_exempt'][0],
		'water_gap':10,
		'shape':'saddle',
		'shape_details':{
			'saddle':{
				'lx':50,
				'ly':50,
				'lz':50,
				'height':12,
				'width':10,
				'binsize':0.9,
				'mono_offset':1.5,
				'lipid':'DPPC',
				},
			'buckle':{
				'lx':40,
				'ly':20,
				'lz':50,
				'height':6,
				'binsize':0.9,
				'mono_offset':1.5,
				'lipid':'DPPC',
				},
			},
		'restraints':{
			'buckle':{'type':'none'},
			'saddle':{'type':'lower_monolayer_extrema','pole_restrain_cutoff':5},
			},
		'top_includes':[
			'martini.ff/martini-v2.2.itp',
			'martini.ff/martini-v2.0-lipids.itp',
			'martini.ff/martini-v2.2-aminoacids.itp',
			'martini.ff/martini-v2.0-ions.itp',
			'martini.ff/martini-v2.0-cholesterol.itp',
			'lipids-tops/martini-v2.0-lipids-zrestrain.itp',
			],
		},
	}
