#!/usr/bin/env tclsh

#---requires tachyon via sudo zypper install tachyon

source sources/codebase/cg_bonds.tcl
source sources/codebase/view_change_render.tcl
color Display Background white
display resize 800 800
display projection orthographic

mol new s01-build-lipidgrid/place-grid-start0.gro
set ind [molinfo top]
mol modstyle 0 $ind Licorice 2.3
set viewpoints(0,$ind,0) { {{0.686704 -0.726339 -0.0295325 0} {0.396589 0.340281 0.852601 0} {-0.609228 -0.597196 0.521731 0} {0 0 0 1}} }
set viewpoints(0,$ind,1) { {{1 0 0 -45.728} {0 1 0 -48.3013} {0 0 1 -15.5706} {0 0 0 1}} }
set viewpoints(0,$ind,2) { {{0.0155263 0 0 0} {0 0.0155263 0 0} {0 0 0.0155263 0} {0 0 0 1}} }
set viewpoints(0,$ind,3) { {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}} }
retrieve_vp 0
render Tachyon scene {
	tachyon -aasamples 12 %s -format TARGA -o %s.tga; 
	convert %s.tga %s.jpeg; 
	rm %s.tga; rm %s;
	}
mol delete 0

mol new s02-build-bilayer/vacuum.gro
cg_bonds -tpr s2-build-bilayer/em-vacuum-steep.tpr
mol modstyle 0 1 Licorice 2
set viewpoints(0,$ind,0) { {{0.686704 -0.726339 -0.0295325 0} {0.396589 0.340281 0.852601 0} {-0.609228 -0.597196 0.521731 0} {0 0 0 1}} }
set viewpoints(0,$ind,1) { {{1 0 0 -45.728} {0 1 0 -48.3013} {0 0 1 -15.5706} {0 0 0 1}} }
set viewpoints(0,$ind,2) { {{0.0155263 0 0 0} {0 0.0155263 0 0} {0 0 0.0155263 0} {0 0 0 1}} }
set viewpoints(0,$ind,3) { {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}} }
retrieve_vp 0
render Tachyon scene2 {
	tachyon -aasamples 12 %s -format TARGA -o %s.tga; 
	convert %s.tga %s.jpeg; 
	rm %s.tga; rm %s;
	}

quit

