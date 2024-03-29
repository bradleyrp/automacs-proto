#!/usr/bin/python

#---SIMULATION RESOLUTION
#---note specify either 'cgmd' or 'aamd' to set the correct parameters
simscale = 'aamd'

#---CONSTRUCTION SETTINGS
construction_settings = {
	'cgmd':{
		#---structures
		'struct_protein':'repo/1ZWW.pdb',
		'top_protein':'repo/snx9/Protein_A.itp',
		'struct_membrane':'repo/membrane-v630/md.part0003.gro',
		'top_membrane':'repo/membrane-v630/system.top',
		#---lattice specs
		'method_lattice':['square','hex'][0],
		'spacer':4,
		'noligos_square':[1,1],
		'noligos_hex':[3,2,3],
		'spacer_within':4.,
		'jitter':[0.,0.],
		'zoffset':0.6,
		#---adhesion method type
		'method_adhesion':['ready','contact','partner'][1],
		#---contact specs
		'struct_lipid':['sources/cgmd-protein-bilayer/orient-PIP2.gro',None][-1],
		'struct_lipid_head_atoms':[1,2,3,4],
		'protein_binding_pocket_residues':[48,91,423,458],
		#---solvent specifications
		'sol_name':'W',
		'negative_ion_name':'CL-','nq':-1,
		'positive_ion_name':'NA+','pq':1,
		'ion_strength':0.150,
		'water_gap_angstroms':5,
		#---bilayer structure characterization
		'director_cgmd':['name PO4','name C4A','name C4B'],
		'selector_cgmd':'name PO4',
		'top_includes':[
			'martini.ff/martini-v2.2.itp',
			'martini.ff/martini-v2.0-lipids.itp',
			'lipids-tops/PIP2.itp',
			'martini.ff/martini-v2.2-aminoacids.itp',
			'martini.ff/martini-v2.0-ions.itp',
			],
		}
	}
	
