# This is an example folder containing the crystal structures of ice Ih and ice III generated using Genice2/2.1.5 program (https://github.com/vitroid/GenIce). 
  We have generated 15 different hydrogen disordered ice structures in order to determine the directional rings dependent upon the proton disordered state.


# USAGE:
1. These options help to count the number of rings formed by at most 10 H-bonds based on H-bond angle pruning criteria.
-> For iceIII
$ nucleation_tracker.py -f iceIII.gro -r 10 -m hbondAngle
Output Files: iceIII_CA_ringsCount.dat

2. These options help to count the number of directional rings formed by at most 10 H-bonds based on H-bond angle pruning criteria.
-> For iceIII
$ nucleation_tracker.py -f iceIII.gro -r 10 -m hbondAngle -d
Output Files: iceIII_DCA_ringsCount.dat
 

3. These options help to count the number of unprunned rings formed by at most 10 H-bonds. 
-> For iceIII
$ nucleation_tracker.py -f iceIII.gro -r 10 -c 1 
Output Files: iceIII_unprunned_ringsCount.dat

4. These options help to count the number of unprunned directional rings formed by at most 10 H-bonds.
-> For iceIII
$ nucleation_tracker.py -f iceIII.gro -r 10 -c 1 -d
Output Files: iceIII_D_unprunned_ringsCount.dat
 

NOTES:
o Unprunned rings mean the rings enumerated without using any prunning criteria/methods (-c 1).
o The 'Fractional directional rings' are calculated from the unprunned rings.

