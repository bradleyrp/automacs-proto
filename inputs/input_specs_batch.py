#!/usr/bin/python -B

#---BATCHMAKER

"""
We create parameter sweeps over simulations by instructing the ``sweep`` variable to sweep over any variable
defined in other specification files set in the ``inputs`` folder. Before we explain how to describe the 
sweep, here is an example workflow for simulating large batches of point mutations from a stock PDB structure. 
In this example, we generate many starting structures and then create a sweep of aamd-protein simulations,
each of which uses a different starting structure.

1. Clone automacs into a directory given by ``batchdir``. If the ``callsign`` variable is a non-empty string then it is attached to all of the downstream runs for easy lookups.
2. Run ``make script protein-homology`` to generate many homology models using a ``sequence.txt`` file in repo or use the mutator to do likewise. If you have a large batch of artisanal protein structures from another source, you can just add them somewhere in the (safe) repo directory and write a list of their pathnames to repo/batch_file_list.txt (referenced below in the sweep).
3. If you use the homology routine, it will write PDB files in subdirectories of s01-homology.
4. Running script-protein-homology will also generate repo/batch_file_list.txt via script-best-models.py, however you may wish to choose the best models yourself and change batch_file_list.txt.
5. Set the sweep variable in ``inputs/input_specs_batch.py`` so that it sweeps over the 'input_filename' variable. 
6. Run ``make batch``
7. This creates a batch directory with subfolders of the form *runN*, numbered over the homology models in ``s01-homology``. The code automatically runs ``make script aamd-protein`` in these folders in order to generate the scripts. 
8. Having automatically generated the code and the necessary script, you must execute the scripts yourself. Note that if you move the data to the cluster you will probably have to run "make rescript" to prepare a cluster script.

The code checks to see if you're doing a batch of mutations by looking at the input_filename variable. If it
looks like a mutation, then it also creates ``batch_names.py`` which gives the user the option to rename
all of the run directories by the name of their mutation. Otherwise, all downstream analysis should refer to
``batch_details.py`` which contains all of the parameters for each run.
   
How do parameter sweeps work?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `sweep` variable is a list of dictionaries, each with "route" and a "values" entries. The route tells 
you what to change, and the values entry tells you how to change it. Routes are lists of keys required to 
re-write any variable in the inputs files. The variable can be buried in a nested dictionary; the top level 
describes variable names defined in th global namespace of the inputs file. In one of the examples below, we 
sweep over the "wbuffer" variable in a protein simulation. The values entry of our sweep hash is normally a 
list of values. The batchmaker will use each of these to create a distinct simulation. If there are multiple 
sweeps, then the batchmaker will make a simulation for each combination of sweep values. For the protein 
mutation example, we have added another option or the values entry. If you use a dictionary containing 
"file_list_file" then the program will use each file in that list as the value for e.g. input_filename.

"""

#---root directory for the batch and universal callsign
batchdir = '../braf-mutants-batch/'
callsign = ''

#---which procedure to execute (same names as "make script")
procedure = ['aamd-protein','cgmd-bilayer'][0]

#---arbitrary parameter sweeps
sweep = {
	'batch_mutation':[
		{'route':['input_filename'],'values':{'file_list_file':'repo/batch_file_list.txt'}},
		],
	'examples':[
		{'route':['protein_construction_settings','aamd','wbuffer'],'values':range(1,1+5)},
		{'route':['complist0'],'values':[[i,1,1] for i in [3,4,5]]},
		],
	}['batch_mutation']

#---overrides for default parameters in auto-generate input-specs files
overrides = {}


