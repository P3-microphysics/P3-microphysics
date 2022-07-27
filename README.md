# Predicted Particle Properties (P3) Microphysics Parameterization

This repository contains the  P3 microphysics scheme and associated components. The package includes the main P3 Fortran module and the lookup tables, all of which are necessary to run with P3 in any model.  It also contains the code necessary to generate the lookup tables, a simple 1D kinematic model (interfaced with P3), and other various things.

## Details on respository contents:

### ./code_p3
This directory contains the main P3 module and the codes necessary to geneate the lookup tables, which are pre-computed and accessed during run-time.  The instructions for generating the lookup tables are found in comment blocks in the code.
- create_p3_lookupTable_1.f90
- create_p3_lookupTable_2.f90
- microphys_p3.f90

### ./interfaces
This subdirectory contains code for interefaces with the GEM and WRF models, for specific model versions, that have been modified to include necessary changes for latest P3 code.
- ./interfaces/gem
- ./interfaces/wrf

### ./kin1d
This subdirectory contains a simple kinematic 1D driving model to perform unit testing of P3.  The README file in that subdirectory describes how to set up and run the model.  A Python plotting script to generate time-height plots from the output is provided.
- Makefile
- README
- ./code
- ./levels
- ./lookup_tables
- ./plots
- ./soundings

### ./lookup_tables
This directory contains the lookup tables (text files), generated from the codes in ./code_p3, necessary to run P3.  During integration, they are read once at the beginning of the time step by s/r P3_INIT and stored in memory.  There are two versions of the main lookup table (1), one for 2-moment-ice and one for 3-moment-ice, and one version of lookup table 2.  Note, the version numbers of the tables are for the tables themselves, and are different than the version number of P3.
- p3_lookupTable_1.dat-v5.4_2momI
- p3_lookupTable_1.dat-v5.4_2momI
- p3_lookupTable_2.dat-v5.3

### ./publications
This subdirectory contains PDF versions of the published scientific artiticles that describe the major developments of the P3 scheme.

### ./playlist
This subdirectory contains links to the main P3 theme songs.
- "Livin la vada loca", Ricky Martin: https://www.youtube.com/watch?v=ltRgb4SJ1uk
- "First we take Manhattan", Leonard Cohen: https://www.youtube.com/watch?v=5rhM1i43NK8


### Further details of the P3 scheme and this repository will be added in the near future.

![P3 mascot](https://user-images.githubusercontent.com/69904757/181296805-14d4de0c-319e-4a28-8a06-5df6d6ac1725.png)
