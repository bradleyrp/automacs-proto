Automatic GROMACS (AMX)
=======================

The "automacs" codes are designed to efficiently, 
reproducibly generate molecular dynamics simulations 
of biophysical systems. The code is commented and 
documentation can be found in the docs folder. This
software is maintained on the 
[github](https://github.com/bradleyrp/automacs)/

Requires
--------

	1. GROMACS version 4 or higher.
	2. Python version 2.5 or higher.
	3. VMD and tachyon for rendering.
	4. The python-Sphinx module for regenerating documenation (development only).

    No installation or compilation required. This software is written almost exclusively 
    in Python, with the addition of some Perl code. We have sought to use a minimum of 
    external packages. The code is designed to run a simulation in-place.

Usage
-----

    Command                          | Outcome
	-------------------------------- | -------------
	./controller                     | Help
	./controller reset               | Clear all simulation data
	./controller make (operation)    | Prepare scripts for "operation"
	./controller make_docs           | Rebuild documentation

Contributing
------------

Ryan Bradley. For other developers, you may view the changelog_specific.md which is updated from time to 
time via a ruby gem called katip. If you have the program, running `katip2.0 changelog_specific.md` will 
update this changelog according to the information from git. A changelog for the project is also located at
changelog.md

