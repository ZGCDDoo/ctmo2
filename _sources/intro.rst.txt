.. _intro:


Introduction
================================

ctmo stands for *Continuous time multi-orbital*.

Sorry the docs are pretty bad for now.

Conventions
----------------

Comments
^^^^^^^^^

We use the convention "$" for the start of a shell command and "#" for commenting shell commands, etc.
    Ex:
        $ cd path/to/thing    # change to the path of thing.



Executable names
^^^^^^^^^^^^^^^^^
There are two main executables: *ctmo* and *ctmo_DCA*. The first executable uses the Cluster Dynamical mean field theory, whereas the second
executable uses Dynamical mean field theory.

Hamiltonian
^^^^^^^^^^^^
The hamiltonian is a density-density type hamiltonian with phonon coupling (**Notice the + sign for the hopping term**):

.. math::

    H  =  \sum_{\rho \rho' \sigma} t_{\rho \rho'} c^\dagger_{\rho\sigma} c_{\rho'\sigma} + \sum_{r, \nu} U_{\nu, \nu} n_{r \nu \uparrow} n_{r \nu \downarrow} + \sum_{r \sigma, \nu < \nu'} U'_{\nu, \nu'} n_{r \nu \sigma} n_{r \nu' \bar{\sigma}}  + 
    \\
    \sum_{r \sigma, \nu < \nu'} \left( U'_{\nu, \nu'} - J^H_{\nu, \nu'} \right) n_{r \nu \sigma} n_{r \nu' \sigma}  + H_{\text{Phonons} }



Important notes for the actual implementation:
    
    - For now, we permit a relatively generic hopping matrix.
    - The value of the interactions are not general yet, only one U, U' and J_H are permitted for now.
    - Only one value (independant of the orbital index) is permitted for the phonon frequency and electron-phonon coupling.