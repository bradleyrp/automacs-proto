CHANGELOG for Automatic GROMACS (AMX)
=====================================

1. **2014.09.23** Major changes to functionality based on user feedback. These include implementation of a make interface which controls controller; removing the automatic restart function from the vacuum backing step in the `amd/bilayer` module and making it into a self-contained tool in `amx/tools`; more careful control of steps and checks in the run control for `amx/bilayer`; and beginning of development of `amx/protein`.