Script started on Fri 28 Nov 2014 11:28:15 AM EST
7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> ls
[0m[01;34mamx[0m  [00mcontroller.py[0m  [01;34mdocs[0m  [01;34minputs[0m  [00mmake[0m  [00mmakefile[0m  [00mREADME.md[0m  [00;32msettings[0m  [01;34msources[0m
7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> more makefile
#---INTERFACE TO CONTROLLER
#------------------------------------------------------------------------------------------------------
-------

#---valid function names from the python script
TARGETS := $(shell perl -n -e '@parts = /^def\s+[a-z,_]+/g; $$\ = "\n"; print for @parts;' \
	controller.py | awk '{print $$2}')

#---collect arguments
RUN_ARGS := $(wordlist 1,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
$(eval $(RUN_ARGS):;@:)

#---hack to always run makefile interfaced with python
scripts=controller.py
$(shell touch $(scripts))
checkfile=.pipeline_up_to_date

#---targets
$(checkfile): $(scripts)
	touch $(checkfile)
	python controller.py ${RUN_ARGS} ${MAKEFLAGS}

#---default and arbitrary make targets
default: $(checkfile)
$(TARGETS): $(checkfile)
[7m--More--(92%)[27m[K
#---git push
gitpush:
	bash amx/script-gitpush.sh ${RUN_ARGS}

7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs>  [Kmore makefilels[Krm -rf repo/[7Pdu -h[3Plscd automacs/l[Ksmv repo-* dev-automacs/ls[K[Kcd automacs/ls[Kdu -hrm -rf repo/ls[Kmore makefile[Kmore makefilels[Krm -rf repo/[7Pdu -h[3Plscd automacs/l[K[Kls
[0m[01;34mamx[0m  [00mcontroller.py[0m  [01;34mdocs[0m  [01;34minputs[0m  [00mmake[0m  [00mmakefile[0m  [00mREADME.md[0m  [00;32msettings[0m  [01;34msources[0m
7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> du -h
244K	./docs/_build/doctrees
8.0K	./docs/_build/html/_static/js
368K	./docs/_build/html/_static/fonts
92K	./docs/_build/html/_static/css
668K	./docs/_build/html/_static
260K	./docs/_build/html/_modules/amx
272K	./docs/_build/html/_modules
12K	./docs/_build/html/_sources
440K	./docs/_build/html/_images
1.5M	./docs/_build/html
1.8M	./docs/_build
1.8M	./docs
8.0K	./.git/info
12K	./.git/refs/remotes/origin
16K	./.git/refs/remotes
8.0K	./.git/refs/heads
4.0K	./.git/refs/tags
32K	./.git/refs
4.0K	./.git/branches
12K	./.git/logs/refs/remotes/origin
16K	./.git/logs/refs/remotes
8.0K	./.git/logs/refs/heads
28K	./.git/logs/refs
36K	./.git/logs
40K	./.git/hooks
4.0K	./.git/objects/info
8.0K	./.git/objects/7d
12K	./.git/objects/e3
8.0K	./.git/objects/37
8.0K	./.git/objects/2d
8.0K	./.git/objects/1e
12K	./.git/objects/cc
16K	./.git/objects/de
24K	./.git/objects/6a
16K	./.git/objects/09
12K	./.git/objects/ce
8.0K	./.git/objects/bd
8.0K	./.git/objects/87
12K	./.git/objects/41
12K	./.git/objects/dc
12K	./.git/objects/b5
8.0K	./.git/objects/df
8.0K	./.git/objects/c9
8.0K	./.git/objects/0c
12K	./.git/objects/c1
12K	./.git/objects/ed
1.5M	./.git/objects/pack
8.0K	./.git/objects/d0
12K	./.git/objects/66
8.0K	./.git/objects/27
1.7M	./.git/objects
1.9M	./.git
24K	./inputs
52K	./sources/codebase
12K	./sources/docs/supplement
464K	./sources/docs
12K	./sources/cgmd-protein-bilayer
8.0K	./sources/aamd-protein-structs
740K	./sources/charmm36.ff
56K	./sources/aamd-bilayer-lipids-structs
72K	./sources/cgmd-bilayer-construct
16K	./sources/general-equil
8.0K	./sources/cgmd-bilayer-lipids-tops
12K	./sources/aamd-protein-homology
16K	./sources/general-sim-restart
12K	./sources/general-sim
20K	./sources/aamd-bilayer-equil
44K	./sources/cgmd-bilayer-equil
388K	./sources/aamd-bilayer-lipids-tops
36K	./sources/aamd-bilayer-construct
160K	./sources/martini.ff
24K	./sources/aamd-protein-equil
40K	./sources/cgmd-protein-construct
12K	./sources/aamd-protein-construct
24K	./sources/cgmd-protein-equil
20K	./sources/cgmd-bilayer-lipids-structs
2.2M	./sources
304K	./amx
6.2M	.
7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> make script aamd-bilayer
touch .pipeline_up_to_date
python controller.py script aamd-bilayer 


    Automacs (AMX) CONTROLLER.

    make script <simulation> : prepare a simulation and associated scripts 
                               see simulation types below for more details
    make script restart      : generates a CPT or GRO restart from an 
                               input-md-in.mdp in repo folder
    make upload              : prepare an upload list and rsync command for 
                               a continuation if the last step was "sim"
    make upload step=<dir>   : prepare an upload list for any step 
                               including the whole directory
    make rescript            : generate continue scripts for a simulation 
                               moved to a cluster  
    make clean               : reset the folder by deleting simulation data
                               note that the repo folder is exempt 

    simulations:
    -----------
    cgmd-bilayer         : coarse-grained bilayer under MARTINI force field
    cgmd-protein         : coarse-grained protein under MARTINI force field 
                           with structure options
    aamd-bilayer         : bilayer under CHARMM36 force field
    aamd-protein         : protein in a water box using CHARMM27 force field
    protein-homology     : single- or muliple-template homology models with 
                           options for batches of mutations
    multiply nx=N ny=M   : replicate a simulation box
    restart              : generic restart which autodetects a CPT or GRO file
    
	Generating singleton script.

	hostname = None
	writing script: script-master-aamd-bilayer
	execute locally with ./script-master-aamd-bilayer
	see the documentation for details

7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> make script aamd-bilayer[K[K[K[K[K[K[K[K[K[K[K[Kcgmd-protein
touch .pipeline_up_to_date
python controller.py script cgmd-protein 


    Automacs (AMX) CONTROLLER.

    make script <simulation> : prepare a simulation and associated scripts 
                               see simulation types below for more details
    make script restart      : generates a CPT or GRO restart from an 
                               input-md-in.mdp in repo folder
    make upload              : prepare an upload list and rsync command for 
                               a continuation if the last step was "sim"
    make upload step=<dir>   : prepare an upload list for any step 
                               including the whole directory
    make rescript            : generate continue scripts for a simulation 
                               moved to a cluster  
    make clean               : reset the folder by deleting simulation data
                               note that the repo folder is exempt 

    simulations:
    -----------
    cgmd-bilayer         : coarse-grained bilayer under MARTINI force field
    cgmd-protein         : coarse-grained protein under MARTINI force field 
                           with structure options
    aamd-bilayer         : bilayer under CHARMM36 force field
    aamd-protein         : protein in a water box using CHARMM27 force field
    protein-homology     : single- or muliple-template homology models with 
                           options for batches of mutations
    multiply nx=N ny=M   : replicate a simulation box
    restart              : generic restart which autodetects a CPT or GRO file
    
	Generating singleton script.

	hostname = None
	writing script: script-master-cgmd-protein
	execute locally with ./script-master-cgmd-protein
	see the documentation for details

7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> make script cgmd-protein[1Pgmd-protein[1Pmd-proteinamd-proteinamd-protein
touch .pipeline_up_to_date
python controller.py script aamd-protein 


    Automacs (AMX) CONTROLLER.

    make script <simulation> : prepare a simulation and associated scripts 
                               see simulation types below for more details
    make script restart      : generates a CPT or GRO restart from an 
                               input-md-in.mdp in repo folder
    make upload              : prepare an upload list and rsync command for 
                               a continuation if the last step was "sim"
    make upload step=<dir>   : prepare an upload list for any step 
                               including the whole directory
    make rescript            : generate continue scripts for a simulation 
                               moved to a cluster  
    make clean               : reset the folder by deleting simulation data
                               note that the repo folder is exempt 

    simulations:
    -----------
    cgmd-bilayer         : coarse-grained bilayer under MARTINI force field
    cgmd-protein         : coarse-grained protein under MARTINI force field 
                           with structure options
    aamd-bilayer         : bilayer under CHARMM36 force field
    aamd-protein         : protein in a water box using CHARMM27 force field
    protein-homology     : single- or muliple-template homology models with 
                           options for batches of mutations
    multiply nx=N ny=M   : replicate a simulation box
    restart              : generic restart which autodetects a CPT or GRO file
    
	Generating singleton script.

	hostname = None
	writing script: script-master-aamd-protein
	execute locally with ./script-master-aamd-protein
	see the documentation for details

7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> grep POSRES_Water[K[K[K[KATER -R
[35m[Kamx/amxsim.py[m[K[36m[K:[m[K			fp.write('#ifdef [01;31m[KPOSRES_WATER[m[K\n')
Binary file amx/amxsim.pyc matches
7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> du -h
244K	./docs/_build/doctrees
8.0K	./docs/_build/html/_static/js
368K	./docs/_build/html/_static/fonts
92K	./docs/_build/html/_static/css
668K	./docs/_build/html/_static
260K	./docs/_build/html/_modules/amx
272K	./docs/_build/html/_modules
12K	./docs/_build/html/_sources
440K	./docs/_build/html/_images
1.5M	./docs/_build/html
1.8M	./docs/_build
1.8M	./docs
16K	./s04-sim
16K	./s03-equil
16K	./s07-sim
8.0K	./.git/info
12K	./.git/refs/remotes/origin
16K	./.git/refs/remotes
8.0K	./.git/refs/heads
4.0K	./.git/refs/tags
32K	./.git/refs
4.0K	./.git/branches
12K	./.git/logs/refs/remotes/origin
16K	./.git/logs/refs/remotes
8.0K	./.git/logs/refs/heads
28K	./.git/logs/refs
36K	./.git/logs
40K	./.git/hooks
4.0K	./.git/objects/info
8.0K	./.git/objects/7d
12K	./.git/objects/e3
8.0K	./.git/objects/37
8.0K	./.git/objects/2d
8.0K	./.git/objects/1e
12K	./.git/objects/cc
16K	./.git/objects/de
24K	./.git/objects/6a
16K	./.git/objects/09
12K	./.git/objects/ce
8.0K	./.git/objects/bd
8.0K	./.git/objects/87
12K	./.git/objects/41
12K	./.git/objects/dc
12K	./.git/objects/b5
8.0K	./.git/objects/df
8.0K	./.git/objects/c9
8.0K	./.git/objects/0c
12K	./.git/objects/c1
12K	./.git/objects/ed
1.5M	./.git/objects/pack
8.0K	./.git/objects/d0
12K	./.git/objects/66
8.0K	./.git/objects/27
1.7M	./.git/objects
1.9M	./.git
20K	./s09-equil
24K	./inputs
52K	./sources/codebase
12K	./sources/docs/supplement
464K	./sources/docs
12K	./sources/cgmd-protein-bilayer
8.0K	./sources/aamd-protein-structs
740K	./sources/charmm36.ff
56K	./sources/aamd-bilayer-lipids-structs
72K	./sources/cgmd-bilayer-construct
16K	./sources/general-equil
8.0K	./sources/cgmd-bilayer-lipids-tops
12K	./sources/aamd-protein-homology
16K	./sources/general-sim-restart
12K	./sources/general-sim
20K	./sources/aamd-bilayer-equil
44K	./sources/cgmd-bilayer-equil
388K	./sources/aamd-bilayer-lipids-tops
36K	./sources/aamd-bilayer-construct
160K	./sources/martini.ff
24K	./sources/aamd-protein-equil
40K	./sources/cgmd-protein-construct
12K	./sources/aamd-protein-construct
24K	./sources/cgmd-protein-equil
20K	./sources/cgmd-bilayer-lipids-structs
2.2M	./sources
304K	./amx
16K	./s06-equil
12K	./s10-sim
6.3M	.
7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> ls
[0m[01;34mamx[0m            [01;34minputs[0m     [01;34ms03-equil[0m  [01;34ms09-equil[0m                   [00;32mscript-master-cgmd-protein[0m
[00mcontroller.py[0m  [00mmake[0m       [01;34ms04-sim[0m    [01;34ms10-sim[0m                     [00;32msettings[0m
[01;34mdocs[0m           [00mmakefile[0m   [01;34ms06-equil[0m  [00;32mscript-master-aamd-bilayer[0m  [01;34msources[0m
[00mgmxpaths.conf[0m  [00mREADME.md[0m  [01;34ms07-sim[0m    [00;32mscript-master-aamd-protein[0m
7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> make clean
touch .pipeline_up_to_date
python controller.py clean 
files to delete:
 ./script-master-aamd-protein
 ./script-master-aamd-bilayer
 ./gmxpaths.conf
 ./script-master-cgmd-protein
 ./s04-sim/settings.sh
 ./s04-sim/script-md-continue
 ./s04-sim/script-sim.pl
 ./s03-equil/script-equilibrate.pl
 ./s03-equil/settings.sh
 ./s07-sim/settings.sh
 ./s07-sim/script-md-continue
 ./s07-sim/script-sim.pl
 ./s09-equil/script-equilibrate.pl
 ./s09-equil/settings.sh
 ./s09-equil/script-md-continue
 ./s06-equil/script-equilibrate.pl
 ./s06-equil/settings.sh
 ./s10-sim/settings.sh
 ./s10-sim/script-sim.pl
directories to delete:
 ./s04-sim/
 ./s03-equil/
 ./s07-sim/
 ./s09-equil/
 ./s06-equil/
 ./s10-sim/
continue? (y/N) y
confirmed? (y/N) y
7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> ls
[0m[01;34mamx[0m  [00mcontroller.py[0m  [01;34mdocs[0m  [01;34minputs[0m  [00mmake[0m  [00mmakefile[0m  [00mREADME.md[0m  [00;32msettings[0m  [01;34msources[0m
7]2;rpb@dark.site:~/compbio/automacs]1;dark.site8rpb@dark:~/compbio/automacs> du -h
244K	./docs/_build/doctrees
8.0K	./docs/_build/html/_static/js
368K	./docs/_build/html/_static/fonts
92K	./docs/_build/html/_static/css
668K	./docs/_build/html/_static
260K	./docs/_build/html/_modules/amx
272K	./docs/_build/html/_modules
12K	./docs/_build/html/_sources
440K	./docs/_build/html/_images
1.5M	./docs/_build/html
1.8M	./docs/_build
1.8M	./docs
8.0K	./.git/info
12K	./.git/refs/remotes/origin
16K	./.git/refs/remotes
8.0K	./.git/refs/heads
4.0K	./.git/refs/tags
32K	./.git/refs
4.0K	./.git/branches
12K	./.git/logs/refs/remotes/origin
16K	./.git/logs/refs/remotes
8.0K	./.git/logs/refs/heads
28K	./.git/logs/refs
36K	./.git/logs
40K	./.git/hooks
4.0K	./.git/objects/info
8.0K	./.git/objects/7d
12K	./.git/objects/e3
8.0K	./.git/objects/37
8.0K	./.git/objects/2d
8.0K	./.git/objects/1e
12K	./.git/objects/cc
16K	./.git/objects/de
24K	./.git/objects/6a
16K	./.git/objects/09
12K	./.git/objects/ce
8.0K	./.git/objects/bd
8.0K	./.git/objects/87
12K	./.git/objects/41
12K	./.git/objects/dc
12K	./.git/objects/b5
8.0K	./.git/objects/df
8.0K	./.git/objects/c9
8.0K	./.git/objects/0c
12K	./.git/objects/c1
12K	./.git/objects/ed
1.5M	./.git/objects/pack
8.0K	./.git/objects/d0
12K	./.git/objects/66
8.0K	./.git/objects/27
1.7M	./.git/objects
1.9M	./.git
24K	./inputs
52K	./sources/codebase
12K	./sources/docs/supplement
464K	./sources/docs
12K	./sources/cgmd-protein-bilayer
8.0K	./sources/aamd-protein-structs
740K	./sources/charmm36.ff
56K	./sources/aamd-bilayer-lipids-structs
72K	./sources/cgmd-bilayer-construct
16K	./sources/general-equil
8.0K	./sources/cgmd-bilayer-lipids-tops
12K	./sources/aamd-protein-homology
16K	./sources/general-sim-restart
12K	./sources/general-sim
20K	./sources/aamd-bilayer-equil
44K	./sources/cgmd-bilayer-equil
388K	./sources/aamd-bilayer-lip