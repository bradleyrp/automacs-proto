Automatic GROMACS (AMX)
=======================

The "automacs" codes are designed to efficiently, 
reproducibly generate molecular dynamics simulations 
of biophysical systems. The code is commented and 
documentation can be found in the docs folder. This
software is maintained on 
[github](https://github.com/bradleyrp/automacs).

Requires
--------

	1. GROMACS version 4 or higher.
	2. Python version 2.5 or higher.
	3. VMD and tachyon for rendering.
	4. The python-Sphinx module for generating documenation.

    No installation or compilation required. This software is written almost exclusively 
    in Python, with the addition of some Perl code. We have sought to use a minimum of 
    external packages. The code is designed to run a simulation in-place.

Usage
-----

Run ``make`` to get started or run ``make docs`` to view the documentation.
