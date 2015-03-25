#!/usr/bin/python

#---COMPOSITIONS
nlipids0 = 600
complist0 = [3,1,1]
lipidnames0 = ['DOPC','DOPS','CHOL']
nlipids1 = 600
complist1 = [3,1,1]
lipidnames1 = ['DOPC','DOPS','CHOL']

#---SIMULATION RESOLUTION
simscale = 'cgmd'

#---CONSTRUCTION SETTINGS

'''
lbuffer : spacing between periodic images
tbuffer : spacing between the ends of lipid tails when making a bilayer from a monolayer
pbc_spacing : distance between monolayers under PBC during vacuum construction step
solvent_thickness : thickness of the water slab in nm, normal to the bilayer (may relax)
lipid_water_buffer : gap between lipids and water during solvate step
ion_strength : molar ion concentration (check the code for the exact procedure)
lipid_spacing : the distance in nm between adjacent lipids during combinator and combinator_row methods
'''

bilayer_construction_settings = {
	'aamd': {
		'lbuffer':0.1,
		'tbuffer':0.1,
		'zbuffer':2.0,
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
		'solvent_thickness':28,
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
		'water_conf':'solvate-water.gro',
		'ion_residue_names':['NA+','CL-','ION'],
		'concentration_calc':['simple','lipid_exempt'][0],
		'top_includes':[
			'martini.ff/martini-v2.2.itp',
			'martini.ff/martini-v2.0-lipids.itp',
			'lipids-tops/PIP2.itp',
			'martini.ff/martini-v2.2-aminoacids.itp',
			'martini.ff/martini-v2.0-ions.itp',
			'martini.ff/martini-v2.0-cholesterol.itp',
			],
		}
	}
	