.. _tutorial:
.. include:: <isotech.txt>

Tutorial
=========

For Mp2 and Graham, use the "scractch" directories. If build problems,
please remove the build directory if it exists, then retry :
    
    $ rm -rf build


Graham (Mp2b): Tutorial 1
---------------------------
1. | Connect to Graham (Mp2b):
   | $ ssh -X "user"@graham.computecanada.ca # where "user" is your compute canada/Mp2 Username
2. Ensure you have done the Pre-Steps described in :ref:`installation`.
3. $ salloc  --time=01:00:00 --ntasks=1 --mem-per-cpu=4000
4. Follow th Graham Procedure in :ref:`installation`.
5. $ module reset 
6. $ module load nixpkgs/16.09  gcc/5.4.0 armadillo boost-mpi
7. $ cd test/Simulations/holstein_DMFT_T0.3
8. $ srun ctmo params1.json

 
Home
---------


Home: Tutorial 1 = Holstein dmft
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you installed on your computer
1. $ cd test/Simulations/holstein_DMFT_T0.3
2. $ bash runCDMFT.sh





