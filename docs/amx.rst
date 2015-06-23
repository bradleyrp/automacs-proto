amx package
===========

*Automacs* is short for "automatic GROMACS" and serves as a catchall for a set of tools designed to
automatically generate molecular dynamics simulations of biophysical lipid bilayers with associated proteins.

Codes are exclusively contained in the ``amx`` module. 
The code itself is designed to coexist with the resulting simulation data.
To get started, run ``make`` to see the current simulation options.
Once you choose a simulation, you can set parameters in configuration files found in ``inputs``.
The ``make`` routine will create an automatic script which can be executed directly or submitted to the queue.

Building coarse-grained bilayers
--------------------------------

The ``lipidgrid`` and ``bilayer`` modules can generate coarse-grained bilayers with a custom composition and size.
This section outlines the general procedure, and by extension, the concept behind the automacs codes.
The basic idea is that we write either Python or Perl codes which perform generic steps such as "make a bilayer" or "equilibrate a simulation".
In general, each step is performed by a distinct Python class or Perl script, and these are executed in sequence by a BASH script.
The procedure for building a coarse-grained bilayer has three parts.

1. Make a two-dimensional grid of lipids in order to form monolayers.
2. Assemble two monolayers into a bilayer and add solvent and counterions.
3. Equilibrate the bilayer system and prepare for a production run.

If you run the following command from the automacs root directory, it will generate the standard bilayer construction codes. ::

	make script cgmd-bilayer

This will generate ``script-master-cgmd-bilayer`` which is written in BASH and is composed of four steps. ::

	#---definitions
	step_monolayer=s01-build-lipidgrid
	step_simulation=s04-sim
	step_build=s02-build-bilayer
	input_files=cgmd-bilayer-equil
	step_equilibration=s03-equil

	#---execute python step
	python -c "execfile('./amx/header');
	try: MonolayerGrids(rootdir='s01-build-lipidgrid');
	except Exception,e:
		print str(e);
		traceback.print_exc();
		print 'fail';"
	go=$(tail s01-build-lipidgrid/log-script-master)
	if [[ "$go" =~ fail$ ]];then exit;fi

	#---execute python step
	python -c "execfile('./amx/header');
	try: Bilayer(rootdir='s02-build-bilayer',previous_dir='s01-build-lipidgrid');
	except Exception,e:
		print str(e);
		traceback.print_exc();
		print 'fail';"
	go=$(tail s02-build-bilayer/log-script-master)
	if [[ "$go" =~ fail$ ]];then exit;fi

	#---execute equilibration step
	cd $step_equilibration
	./script-equilibrate.pl
	if [[ $(tail log-script-master) =~ fail$ ]];then exit;fi
	cd ..

	#---execute simulation (continuation) step
	cd $step_simulation
	./script-sim.pl
	if [[ $(tail log-script-master) =~ fail$ ]];then exit;fi
	cd ..


The first two steps instantiate the MonolayerGrids and Bilayer classes.
These classes will access settings found in the ``inputs/input-specs-bilayer.dat`` file in order to design the resulting bilayer.
After constructing the bilayer using the Python codes, we perform equilibration and simulation using standard Perl scripts, which are less flexible but easier to implement on clusters.

The Python classes each create a new folder which represents a discrete step in the construction process.
For example, making a bilayer happens in two steps: first we build two monolayers and then we adhere them into a bilayer.
After that, we use Perl scripts to run equilibration and production simulations in two extra folders. 
The folders are always named by their step number so everything stays in order.

This workflow allows each class instance to operate on a self-contained folder in order to assemble and 
simulate the system of interest. This minimizes the amount and complexity of the code, since we repeat many 
steps (e.g. regenerate a different set of counterions) during the course of running many simulations. 

By convention, the input specifications are stored in the ``inputs`` folder while chemical information, such 
as lipid configurations and topologies are stored in the ``sources`` folder. All of these settings can be 
changed according to the class specifications given below.

amx.tools module
----------------

.. automodule:: amx.tools
    :members:
    :undoc-members:
    :show-inheritance:

