amx Module
===========

All automacs codes are contained in this module. *Automacs* is short for *automatic GROMACS*, and serves as a catchall for a set of tools designed to
automatically generate molecular dynamics simulations of biophysical lipid bilayers with associated proteins.
Bleeding edge features include batch simulation of mutated proteins, more articulated bilayer geometries, and 
more advanced configurations of protein-membrane assemblies.

.. _sec-concept:

Concept
-------

Automacs is designed to eliminate or reduce the amount of repetitive labor required to create a biophysical simulation
in GROMACS. For that reason, we have attempted reduce most human-readable parameters to a small collection of 
input files while automating all of the repetitive steps which are common to many simulations. By focusing on a collection
of tuneable parameters, we hope to make it possible to generate large simulation datasets with a minimum of description
and so-called "manual" labor which invites mistakes and wastes time.

The user controls the software from the ``makefile`` which orchestrates a small collection of Python, Perl, and Bash codes. 
If you don't like reading this documentation, you can simply run ``make`` to see which options are available. In general,
you can start a simulation with the following procedure.

1. Choose a scientific question that you can answer with one of our procedures.
2. Edit one of the python dictionaries found in ``inputs/input_specs_PROCEDURE.py``.
3. Run ``make script PROCEDURE`` which will create an executable called ``script-PROCEDURE`` which you can execute from the terminal.

Tunable parameters can be found in three different forms.

1. The files found in the inputs directory contain nested python dictionaries which set all of the necessary geometries, compositions, and simulation parameters.
2. The settings file in the root directory contains paths, software version numbers, hostnames, cluster execution scripts, processor settings, and other platform-specific details necessary to make sure the software runs correctly on your machine.
3. The ``sources`` folder holds some input files, structures, and simulation scripts.
4. All of the other parameters are written into the codes executed from ``controller.py``. These parameters are only intended for advanced users.

Each procedure follows the same basic workflow.

1. The ``makefile`` calls ``controller.py`` which looks at the appropriate ``inputs/input_specs_PROCEDURE.py`` file to set up the right simulation scripts in the correct order.
2. The automatic script (usually named ``script-PROCEDURE``) consists of a list of discrete steps, each of which corresponds to a distinct directory found in the root automacs directory (sometimes these are generated at run-time). These folders are numbered for readability (e.g. ``s03-equil``).

	a. Python steps are written directly in Bash and start by executing ``amx/header`` which loads our codes. Each step is executed by creating a class object. We usually send along the desired step directory and sometimes the name of the previous directory. All of the other parameters come from the ``inputs`` folder.
	b. Perl steps, preferred for compatibility with clusters, are executed directly from the script, and usually look to a customized ``settings.sh`` file to get the right pathnames and settings. We use Perl to execute equilibration and production run steps. Bash generally passes the correct variables to the downstream codes.
	
In the course of executing these steps, the construction and simulation settings are retrieved from the sources listed above, otherwise they are written directly into the code. We try to use as few general settings as possible; as a result, the inputs files can be fairly descriptive. This workflow allows each class instance to operate on a self-contained folder in order to assemble and simulate the system of interest. This minimizes the amount and complexity of the code, since we repeat many steps (e.g. regenerate the counterions, equilibrate a system) during the course of running many simulations. 

Bilayers
--------

Coarse-grained bilayers
^^^^^^^^^^^^^^^^^^^^^^^

The procedure for building a coarse-grained bilayer has three parts.

1. Make a two-dimensional grid of lipids in order to form monolayers.
2. Assemble two monolayers into a bilayer and add solvent and counterions.
3. Equilibrate the bilayer system and prepare for a production run.

These steps are performed by the ``lipidgrid`` and ``bilayer`` modules which prodice the first two steps in the bilayer construction procedure, normally stored in folders titled ``s01-build-lipidgrid`` and ``s02-build-bilayer``. The resulting structures are simulated by general-purpose Perl codes for equilibration and production runs which generate the ``s03-equil`` and ``s04-sim`` folders.

Lipid grids
"""""""""""

.. automodule:: amx.lipidgrid
    :members:
    :undoc-members:
    :show-inheritance:

Bilayers
""""""""

.. automodule:: amx.bilayer
    :members:
    :undoc-members:
    :show-inheritance:

Homology modeling
-----------------

Rudimentary homology modeling is implemented via `MODELLER <https://salilab.org/modeller/>`_. 
The code has two main uses. You can create a model from one or more templates, or you can make a large 
batch of mutations from either a stock PDB structure or a custom structure that you drop in ``repo``.
The latter can be combined with the :ref:`sec-batchmaker` to run large numbers of protein simulations.

Here is an example procedure in which the user creates a homology model of a known structure, adds a point mutation, and then runs a protein-in-water simulation.

1. Choose a structure. Locate the tuple in the *mutator, template* section of ``homology_construction_settings`` in the ``inputs/input_specs_homology.py``. This tuple has two items: a the "name" of the protein and its chain. The chain is always a single letter which you must specify.  You can name your protein in one of two ways.
	a. If you provide a PDB code (always a 4-character string), the program will download that PDB for you automatically.
	b. If you have a custom structure, add it to the ``repo`` folder, and use the file prefix (i.e. without ".pdb" suffix) as its name.
2. Enter your mutations in *mutator, mutations* which is a list of tuples. Each tuple contains three items: the original residue (we check that this matches), the position, and the mutated residue. We use single-letter amino acid codes.
3. MODELLER will make many models so that you can pick the best one. We recommend setting ``nmodels`` in *mutator* to at least 20.
4. Run ``make script protein-homology``. This runs MODELLER and creates the starting structure.
5. The file ``repo/batch_file_list.txt`` contains the paths to the top model for your mutation(s). Choose one carefully. 
6. Replace ``input_filename`` in ``inputs/input_specs_protein.py`` with your chosen path. This will tell the simulator where to find the starting structure.
7. Procede with an atomistic protein simulation as per the instructions in :ref:`sec-aamd-protein` (see step 3, since you've already set ``input_filename``). 

.. automodule:: amx.proteinhomology
    :members:
    :undoc-members:
    :show-inheritance:

.. _sec-aamd-protein:

Atomistic protein simulation
----------------------------

Creating a simulation of an atomistic protein in water is one of the most common modeling tasks. All of the tunable parameters for this procedure are located in ``inputs/input_specs_protein.py``, however they are set to industry-standard defaults and require little intervention. The user must only supply a protein structure to get started. Here is a simple use case.

1. Choose a structure. Save it to the ``repo`` folder (create this folder if it doesn't exist already).
2. Set the ``input_filename`` variable in ``inputs/input_specs_protein.py`` to the path for your structure. This path is relative to the automacs root directory so you should include the folder name (e.g. ``repo/my_structure.pdb``). The ``repo`` directory is safe from deletion when you reset the code via ``make clean``.
3. Run ``make script aamd-protein`` to prepare the simulator scripts.
4. Run ``./script-aamd-protein`` to run the simulation.

This will produce a standard atomistic simulation. To increase the duration of the simulation, you can run ``./script-md-continue`` in the ``sN-sim`` folder, which contains the production run trajectories. You can continue the simulation as many times as you want.

.. automodule:: amx.protein
    :members:
    :undoc-members:
    :show-inheritance:

.. _sec-batchmaker:

Batchmaker
----------

Automacs is self-replicating. All of the parameters defined in the ``inputs`` folder can be swept over an 
arbitrary number of parameters in combinations. This can be done by setting the ``sweep`` variable in the 
``inputs/input_specs_batch.py`` file.

.. automodule:: inputs.input_specs_batch
    :members:
    :undoc-members:
    :show-inheritance:

Equilibration
-------------

For compatibility purposes, both equilibration and production simulations are executed via Perl scripts.

General codes
-------------

.. automodule:: amx.amxsim
    :members:
    :undoc-members:
    :show-inheritance:

Tools
-----

.. automodule:: amx.tools
    :members:
    :undoc-members:
    :show-inheritance:

