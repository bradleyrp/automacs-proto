#!/bin/bash

#---AUTOMACS (automatic GROMACS) run settings
#---from sources/general-equil/settings.sh

#---equilibration steps (nvt-short, nvt, npt)
STEPSTRING="nvt-short nvt npt"

#---input files
SOURCEDIR=
MDPLOOK=
SPECSFILE=
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
cp $SOURCEDIR/$ITP ./
cp $SOURCEDIR/composition.dat ./
cp -r $SOURCEDIR/*.ff ./
cp -r $SOURCEDIR/lipids-tops ./

#---copy standard input files for minimization (deprecated via input_specs_mdp.py)
breadcrumb=$(pwd);cd ../
python -c "execfile('amx/header');
try: 
	sim = AMXSimulation();
	sim.params = {};
	execfile('"$SPECSFILE"',sim.params)
	sim.write_mdp(mdp_route='"$MDPLOOK"',mdp_section='mdp-equil',rootdir='"$breadcrumb"')
except Exception,e:
	print str(e);
	traceback.print_exc();
	print 'fail';"
cd $breadcrumb
