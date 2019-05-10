.. _tutorial:
.. include:: <isotech.txt>

Tutorial
=========

For Mp2 and Graham, use the "scractch" directories. If build problems,
please remove the build directory if it exists, then retry :
    
    $ rm -rf build


Graham (Mp2b): Tutorial 1
---------------------------


First steps
^^^^^^^^^^^^^^^^^^^^^^^^
1. | Connect to Graham (Mp2b):
   | $ ssh -X "user"@graham.computecanada.ca # where "user" is your compute canada/Mp2 Username
2. Ensure you have done the Pre-Steps described in :ref:`installation`.
3. $ salloc  --time=01:00:00 --ntasks=1 --mem-per-cpu=4000 --account=def-tremblay
4. $ module reset 
5. $ module load nixpkgs/16.09  gcc/5.4.0 armadillo boost-mpi
6. $ cd ~/scratch
7. $ cp -r ~/projects/def-tremblay/tutorial_ctmo ./
8. $ srun ctmo params1.json 

Here are the next steps.

Next steps
^^^^^^^^^^^^^^^
This will produce many files. Look at them with ls:

The important output files are greenUp.dat, hybUp.dat, selUp.dat. Note that the files are renamed corresponding to
the current iteration.

Let us plot the green fonction.

    $ gnuplot
    $ plot 'greenUp.dat' u 1:3 w lp ps 2



