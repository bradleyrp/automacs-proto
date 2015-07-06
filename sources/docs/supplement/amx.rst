amx package
===========

"Automacs" is short for automatic GROMACS, and serves as a catchall for a set of tools designed to
automatically generate molecular dynamics simulations of biophysical lipid bilayers with associated proteins.

Currently, the ``amx`` package contains only routines to build MARTINI coarse-grained bilayers. 

Constructing CGMD bilayers
--------------------------

The procedure for building a coarse-grained bilayer has three parts.

1. Make a two-dimensional grid of lipids in order to form monolayers.
2. Assemble two monolayers into a bilayer and add solvent and counterions.
3. Equilibrate the bilayer system and prepare for a production run.

Assembly of the monolayer and bilayer are performed by the ``lipidgrid`` and ``bilayer`` modules, 
respectively. These modules are described in detail below. These codes are executed by the python script 
titled ``script-amx-cgmd-bilayer``. This script contains the following key commands. ::

	from amx.lipidgrid import MonolayerGrids
	from amx.bilayer import Bilayer
	from amx.tools import *

	MonolayerGrids(rootdir='s1-build-lipidgrid')
	Bilayer(rootdir='s2-build-bilayer',previous_dir='s1-build-lipidgrid')

After importing the necessary classes from the ``amx`` module, the script creates a ``MonolayerGrids`` 
instance followed by a ``Bilayer`` instance. These class instances will generate simulation data in new 
subfolders, which are preferably labeled by step. That is ``s2-build-bilayer`` refers to the second step in 
which we build the bilayer from the monolayers generated in the first. 

This workflow allows each class instance to operate on a self-contained folder in order to assemble and 
simulate the system of interest. This minimizes the amount and complexity of the code, since we repeat many 
steps (e.g. regenerate a different set of counterions) during the course of running many simulations. 

By convention, the input specifications are stored in the ``inputs`` folder while chemical information, such 
as lipid configurations and topologies are stored in the ``sources`` folder. All of these settings can be 
changed according to the class specifications given below.

Equilibration
-------------

For compatibility purposes, both equilibration and production simulations are executed via Perl scripts.

General codes
-------------

.. automodule:: amx.amxsim
    :members:
    :undoc-members:
    :show-inheritance:


Lipid grids
-----------

.. automodule:: amx.lipidgrid
    :members:
    :undoc-members:
    :show-inheritance:

Bilayers
--------

.. automodule:: amx.bilayer
    :members:
    :undoc-members:
    :show-inheritance:

Tools
-----

.. automodule:: amx.tools
    :members:
    :undoc-members:
    :show-inheritance:

