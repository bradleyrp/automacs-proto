#!/bin/bash

#---AUTOMACS (automatic GROMACS) run settings
#---from sources/general-sim/settings.sh

#---input files
SOURCEDIR=../s2-equil
EXTENDTIME="50000"
STOPHOURS="-1"

if ls *cpt 1> /dev/null 2>&1; then
echo "settings.sh found a cpt file so skipping copy procedure"
else

#---locate the latest checkpoint from the equilibrate folder
THISDIR=$(pwd)
cd $SOURCEDIR
PRUN=0
for file in md.part*.cpt
do
	if [ $(echo ${file:7:4} | sed 's/^0*//') -gt $PRUN ]; then 
		PRUN=$(echo ${file:7:4} | sed 's/^0*//')
	fi
done
cd $THISDIR

#---copy checkpoint and input files
cp $SOURCEDIR/$(printf md.part%04d.tpr $PRUN) ./
cp $SOURCEDIR/$(printf md.part%04d.cpt $PRUN) ./

fi
