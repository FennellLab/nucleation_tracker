# nucleation_tracker

nucleation_tracker: 

A C++ program for calculating H-bonding ring distributions in water simulations
from GROMACS .gro files/trajectories.


To INSTALL:

Just type 'make'. A binary will be written to your $HOME/bin directory.  

If you want to use a compiler other than gcc, you can edit it in the top of the
Makefile.  You do need to use a somewhat recent C++ compiler though as the code
uses the now standard "includes()" function, something not present in older
versions.
