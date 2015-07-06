#!/bin/bash

#---AUTOMACS (automatic GROMACS) run settings
#---from sources/flexible-equil/settings.sh

#---input files
SOURCEDIR=
STRUCT=system.gro
TOP=system.top
ITP=*.itp
GRP=system-groups.ndx

#---path to input repository
AMXPATH=../sources

#---copy input files from other steps
cp $SOURCEDIR/$STRUCT ./system-input.gro
cp $SOURCEDIR/$TOP ./system.top
cp $SOURCEDIR/$GRP ./system-groups.ndx
cp $SOURCEDIR/$ITfdasP ./
cp -r $SOURCEDIR/*.ff ./
cp -r $SOURCEDIR/lipids-tops ./

#---needs copied from general-equil
#---copy standard input files for minimization (deprecated via input-specs-mdp.dat)
#---execute python step
python -c "execfile('../amx/header');
try: 
	sim = AMXSimulation();
	sim.params = {};
	execfile('inputs/input-specs-bilayer.dat',sim.params)
	sim.write_mdp(mdp_section='mdp-equil',rootdir='"$SOURCEDIR"')
except Exception,e:
	print str(e);
	traceback.print_exc();
	print 'fail';"

