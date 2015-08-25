#!/usr/bin/python -i

"""
CUSTOM AUTOMACS SETTINGS 
system-specific execution settings

Note that this file is a good example of what ~/.automacs.py should look like. 

This file adds definitions to pre-defined automacs globals which are initially set in the settings file in 
the automacs root. If you don't add any settings below, automacs will call GROMACS using default
settings. For example, it will call "mdrun" directly. If you rely on the "modules" system to load your 
GROMACS paths, or you wish to change the number of threads, set up a cluster script, etc, then you must
tell automacs below.

Custom settings are used only if the system's hostname matches a key in `valid_hostnames`. The hostname is 
checked either by echoing the bash variable of the same name or by using the socket module in python as a 
failsafe. The key can either be the explicit hostname or a regex that matches it (e.g. you can use "cluster" 
even if the hostname is "cluster.login"). 

If you set valid_hostnames[hostname_key] equal to `None` then you can define one set of custom settings 
for your system. If you want more than one setting in the case of multiple architectures, different PBS 
queues, and so on, then you should set valid_hostnames[hostname_key] to some other label. This label will 
be the default for that system, and you can change it to submit to e.g. different queues. Once you choose a 
label, you can set the other variables using an underscore-concatenated key equal to 
`<hostname_key>_<label>`. For example, if we want to use the xeons on cluster.login, we might set 
`valid_hostnames['cluster'] = 'xeons'` and substitute `cluster_xeons` as the system_name in all of the 
examples below.

There are two dictionaries which override the standard GROMACS commands: gmx_overrides and gmx_suffixes.

If you are running on a system that requires `mpirun` to run the GROMACS `mdrun` utility, then you can add an 
entry to the gmx_overrides dictionary using the common name for that gromacs utility. Each entry should be a 
dictionary whos keys are standard GROMACS commands. This also lets you set the number of threads or 
processors. For example, if we wish to run `mdrun` with four processors using mpirun, then we might use: 
`gmx_overrides[system_name] = {'mdrun':'mpirun -np 4 mdrun'}`. This would also allow us to force the use the 
GPU if one is available using `gmx_overrides[system_name] = {'mdrun':'mdrun -nb cpu'}`.

If all of your GROMACS executables were generated with a suffix e.g. `mdrun_mpi_d` then you can add an entry 
to gmx_suffixes to append that string to all command-line executables (e.g. 
`gmx_suffixes[system_name] = '_mpi_d'`).

For running on clusters, we add a sub-dictionary to default_proc_specs. The resulting 
default_proc_specs[system_name] dictionary can have the following entries.

1. nodes : the number of nodes
2. ppn : the number of processors per node
3. walltime : number of hours of walltime to request
4. module : a newline-separated string of all of the `modules` commands to run to load the right GROMACS paths

Note that you may use `NPROCS` in gmx_overrides if you want automacs to compute the correct number of threads 
or processors to use in any `mdrun` commands if you're running the simulation on a cluster. The settings 
listed above will only apply if you also define a variable called `cluster_header_system_name` which must 
contain the necessary lines. Automacs will use the default_proc_specs dictionary to override these settings, 
but the lines have to be present in the cluster header. Here is an example.

    cluster_header_kraken_xeons = \
    \"\"\"#!/bin/bash
    #PBS -l nodes=1:ppn=16,walltime=96:00:00
    #PBS -j eo 
    #PBS -q opterons
    #PBS -N gmxjob
    echo "Job started on `hostname` at `date`"
    cd $PBS_O_WORKDIR
	\"\"\"

To reiterate: you have to have a nodes/ppn/walltime line in order to use default_proc_specs. Module loading 
commands are appended to the cluster header.

Finally, if you are running GROMACS 5 or later, automacs will try to run `gmx` and if it succeeds without 
error, then it assumes that a 5.0 series GROMACS is available, in which case it uses the slightly different 
utility names, though the remainder of the system works the same.
"""


#---a local machine with a GPU and only one ar
valid_hostnames['dark'] = None
gmx_overrides['dark'] = {'mdrun':'mdrun -nb cpu'}

#---a cluster with a PBS header and a particular architecture
valid_hostnames['compbio'] = 'opterons'
default_proc_specs['compbio_opterons'] = {
	'nodes':1,'ppn':16,'walltime':48,'module':'module load gromacs-gcc','scratch':False}
cluster_header_compbio_opterons = \
"""#!/bin/bash
#PBS -l nodes=1:ppn=16,walltime=96:00:00
#PBS -j eo 
#PBS -q opterons
#PBS -N gmxjob
echo "Job started on `hostname` at `date`"
cd $PBS_O_WORKDIR
"""

