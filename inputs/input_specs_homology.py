#!/usr/bin/python

"""
HOMOLOGY MODELING routines based on MODELLER.

This code wraps MODELLER so that you can predict protein structures or make mutations. These codes are capable 
of creating protein structures from one or many templates. The 'mutator' function also allows you to make 
large batches of point mutations to proteins in the PDB.

Templates are specified by tuples that contain the PDB code and the correct protein chain (usually a letter).
Targets are specified in one of several different formats. We use a dictionary for clarity.

1. A raw sequence input containing a tuple with the protein name and sequence.
2. A text file (in the safe "repo" directory) which provides a list of sequences.
3. A text file with a list of RNA sequences.
4. A verbatim template from the PDB (which might be missing part of the sequences.

Here are some examples of these four types:

'target':{'raw':('mdia2NB','MERHQPRLHHPAQGSAAGTPYPSSASLRGCRESKMPRRKGPQHPPPPSGPEEPGEKRPKFHLNIRT')},
'target':{'textfile':'repo/sequences.txt'},
'target':{'textfile_rna':'repo/sequences_rna.txt'},
'template':[('1LP1','A')],

Note that target_name must always be a list in the 'single' procedure. If target_name is not defined, the 
'single' procedure will make a batch of models.
"""

#---choose the modeller procedure detailed in homology_construciton settings
procedure = [
	'single',
	'multi',
	'mutator',
	][2]

homology_construction_settings = {

	#---predict protein structure from a single PDB template
	'single':{
		'target':{'raw':('plcdPH',
			'GLQDDPDLQALLKGSQLLKVKSSSWRRERFYKLQEDCKTIWQESRKVMRSPESQLFSIED'+
			'IQEVRMGHRTEGLEKFARDIPEDRCFSIVFKDQRNTLDLIAPSPADAQHWVQGLRKIIH')},
		'regex_subs':[('tttatttagagcctg','tttattagcctg')],
		'template':[('1MAI','A')],
		'n_models':10,
		},

	#---predict protein structure from multiple PDB templates
	'multi':{
		'target':{'textfile':'repo/sequences.txt'},
		'template':[
			('3O4X','A'),
			('2RCA','B'),
			],
		'n_models':2,
		},

	#---mutate a protein from the PDB
	'mutator':{
		'example_from_custom':{
			'n_models':4,
			'template':[
				('braf_inactive','B'),
				],
			'mutations':[
				('B',468+1,'J'),
				][:],			
			},
		'example_from_pdb':{
			'n_models':2,
			'template':('2GS7','A'),
			'mutations':[
				('K',690,'R'),
				('E',710,'R'),
				('E',710,'A'),
				],
			},
		}['example_from_pdb'],
	}

