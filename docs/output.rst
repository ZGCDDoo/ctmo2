.. _output:
.. include:: <isotech.txt>


Conventions
=============
The conventions for the cluster sites follow the cartesian plan, starting by the site 0.
The sites are numbered in the order of increasing x, then increasing y, then increasing z.

The same applies to the definition of the K points in DCA.

Important Output files
=======================
The important output files are the following:

    - outPutConvention.dat:
        The outPut of the green, hyb and self files regarding the cluster sites.
        For example, if the outPutConvention file has the following form:
        (0,0)
        (0,1)
        (0,3)
        
        We thus have N Independant green fuctions (NIndep "Sites pair").
        Then the green function has the first column as the frequencies, the second and third column
        as the real and imaginary part of G_00 respectively and so on.
        If there are multiple orbitals, then we repeat this process for each orbital.
        Ex. Say 2 orbitals, o1 and o2, then we have the previous definition for (o1, o1), then (o1, o2), then (o2, o2)
        for a total of (3*3*2 + 1) columns, i.e 1 + [NIndep*NOrb*(NOrb+1)/2 ]*2 Columns.
        

    - greenUp${ITER}.dat:
        The green function.

    - selfUp${ITER}.dat:
        The self-energy function.


    - hybUp${ITER}.dat:
        The hybridization function.

    - Obs.json:
        All the observables computed. For each of these observables, we also have
        a .dat file associated. For example, docc.dat, sign.dat, docc0_0.dat (double occupancy of site0, orbital0).
        Etc.


Averaging the results
=====================

Simple as using a python module:

1. module load python/3.7 scipy-stack
2. pip install --user statsfiles
3. cd DIRECTORY_OF_SIMULAION
4. Replace ${ITERSTART} in the following line with the iteration at which
    to start the averging.
5. python -m statsfiles ${ITERSTART} --cdh


Sometimes, you may encounter problems if the program ctmo has been terminated while
it would have been writing files, such as the green function.
Therefore, if you get an error of numpy shape, then remove the last iteraion of greenUp, selfUp and hybUp
files and this should resolve the problem.

