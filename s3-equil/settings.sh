#!/bin/bash

#---AUTOMACS (automatic GROMACS) run settings

#---type
SIMTYPE="aamd-bilayer"

#---equilibration steps (nvt-short, nvt, npt)
STEPSTRING="nvt-short nvt npt"

#---input files
SOURCEDIR=../s2-build-bilayer
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
cp -r $SOURCEDIR/martini.ff ./
cp -r $SOURCEDIR/lipids-tops ./

#---copy standard input files for minimization
cp $AMXPATH/aamd-bilayer-equil/* ./
