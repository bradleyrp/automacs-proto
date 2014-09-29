#!/bin/bash

#---AUTOMACS (automatic GROMACS) run settings

#---type
SIMTYPE="aamd-protein"

#---equilibration steps (nvt-short, nvt, npt)
STEPSTRING="nvt npt"

#---input files
SOURCEDIR=../s1-build-protein-water
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

#---copy standard input files for minimization
cp $AMXPATH/aamd-protein-equil/* ./
