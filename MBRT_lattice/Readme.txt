The code and data related to minibeam lattice are available here: https://drive.google.com/drive/folders/1PRfH6sstmxG0weiGKZe4r4nsL-4MrB4Q?usp=sharing
Few key points:
1. Inside the 'matRad-master-collimator' folder: there is will be 'genDij_9306087_10212024.m' and 'genDij_2286842_10212024.m' file. You can use this to generate Dij's for minibeams. Specifically, in this file, 
(a) set grid dimension in line 87 using the variable 'resolution'. Note that the default resolution is 3x3x3 mm. 
(b) Set center to center distance between slits of multi slit collimator on line 65 using variable 'cs_coll'. 
(c) Set orientation of collimator wrt to beam source in line 66 using variable 'cs_angle'. This refers to the angle (in degrees). Current 'genDij...' can only handle values of 0 or 90 degrees here. 
(d) Set the angle using the variable 'pln.propStf.gantryAngles' on line 37.

2. Please refer to the 'test_MBL_DO_9306087_011224.m' and 'test_MBL_DO_2286842_011224.m' files to run the algorithm for the two test cases. 
(a) Starting on line 26 in both files, we define the lattice peaks, i.e., the x,y,z coordinates of lattice peaks are defined using 'lattice_x', 'lattice_y', and 'lattice_z' variables. 'ddr' is the diamter of the peak.
(b) Please change the 'folder' variable as well as path to the files at the beginning of the code appropriately.
(c) On line 4, 'method' variable is used to set the method that is used. Please refer to the paper for the details about the methods.