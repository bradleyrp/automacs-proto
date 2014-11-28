#!/bin/bash

#---AUTOMACS (automatic GROMACS) run settings

#---equilibration steps (nvt-short, nvt, npt)
STEPSTRING="nvt-short nvt npt"

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
cp $SOURCEDIR/$ITP ./
cp $SOURCEDIR/composition.dat ./
cp -r $SOURCEDIR/*.ff ./
cp -r $SOURCEDIR/lipids-tops ./

#---copy standard input files for minimization
cp $AMXPATH/cgmd-bilayer-equil/* ./
