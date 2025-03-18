# This is an example folder containing a multi-frame trajectory file (liqWat_T4I.gro) of pure water MD simulation at 298.15 K temperature using TIP4P/Ice water model.


# USAGE:
1. These options help to count the number of rings formed by at most 10 H-bonds based on H-bond angle pruning criteria.
$ nucleation_tracker.py -f liqWat_T4I.gro -r 10 -m hbondAngle
Output Files: liqWat_T4I_CA_ringsCount.dat

2. These options help to count the number of directional rings formed by at most 10 H-bonds based on H-bond angle pruning criteria.
$ nucleation_tracker.py -f liqWat_T4I.gro -r 10 -m hbondAngle -d
Output Files: liqWat_T4I_DCA_ringsCount.dat
 

3. These options help to count the number of unprunned rings formed by at most 10 H-bonds. 
$ nucleation_tracker.py -f liqWat_T4I.gro -r 10 -c 1
Output Files: liqWat_T4I_unprunned_ringsCount.dat

4. These options help to count the number of unprunned directional rings formed by at most 10 H-bonds.
$ nucleation_tracker.py -f liqWat_T4I.gro -r 10 -c 1 -d
Output Files: liqWat_T4I_D_unprunned_ringsCount.dat
   

NOTES:
o Unprunned rings mean the rings enumerated without using any prunning criteria/methods (-c 1).
o The 'Fractional directional rings' are calculated from the unprunned rings.
o The RSF order parameter and tetrahedral order parameter are found out using the options '-s' and '-t' respectively. The RSF is output in the main 'ringsCount.dat' file, and a separate folder is created to output the tetrahedral order parameter along with pdb file with all water tetrahedrality.
o The POVRay files are generated using '-p' option with a commandline to generate the POVRay images.
o The ring location '.xyz' file is generated using '-x' option, represented by the center of geometry of the polygons. 

