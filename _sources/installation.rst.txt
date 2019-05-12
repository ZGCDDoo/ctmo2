.. _installation:

Installation
================================


**Note :**
If build problems,
please remove the build directory if it exists, then retry :
    
    $ rm -rf build

Dependencies
--------------

Mandatory Dependencies
^^^^^^^^^^^^^^^^^^^^^^^
1. Armadillo
2. boost (mpi, serialization, filesystem, system)

Optional Dependencies
^^^^^^^^^^^^^^^^^^^^^^
1. snappy
2. Postgresql


Pre-Steps
----------
1. Make sure you have a "bin" directory in your home folder
2. Append the bin folder to your path. Add the following line to your ~/.bashrc:  export PATH="$PATH:~/bin"
3. $ source ~/.bashrc

Linux (Ubuntu 16.04)
----------------------
This installation procedure should work for many recent Linux flavors. For the following
we present the instructions specific for Ubuntu or derivatives.

1. Install the Dependencies
    $ sudo apt-get install libarmadillo-dev libboost-all-dev cmake
2. | $ mkdir build && cd build && cmake -DBUILD_TESTS=OFF .. && make -j NUMBER_OF_CORES
   | # replace NUMBER_OF_CORE by say = 4
3. sudo make install


Mac
-----
This has been tested on macOS 10.13.6.

1. Install the Dependencies (with Homebrew : https://brew.sh/)
    $ brew install armadillo
    $ brew install boost
    $ brew install boost-mpi
    $ brew install postgresql 
    $ brew install snappy
2. | $ mkdir build && cd build && cmake -DBUILD_MAC=ON -DBUILD_HOME=OFF .. && make -j NUMBER_OF_CORES
   | # replace NUMBER_OF_CORE by say = 4
3. sudo make install



Graham, Ceder, mp2b, ms2b
--------------------------
1. $ module reset 
2. $ module load nixpkgs/16.09  gcc/5.4.0 armadillo boost-mpi snappy
3. | $ mkdir build && cd build
   | CXX=mpic++ cmake -DBUILD_GRAHAM=ON -DBUILD_MPI=ON  -DBUILD_HOME=OFF .. 
   | make
4. Copy the executaables to a know location on your path.
   | For example: $ mkdir ~/bin && cp -i ../src/ctmo* ~/bin/
   | Add the ~/bin to your path : 
   | in ~/.bashrc add the following line at the end:
   export PATH="$PATH:~/bin"



Beluga
-------
1. $ module reset 
2. $ module load nixpkgs/16.09 gcc/7.3.0 cmake armadillo/7.950.1 snappy boost-mpi
3. | $ mkdir build && cd build
   | CXX=mpic++ cmake -DBUILD_GRAHAM=ON -DBUILD_MPI=ON  -DBUILD_HOME=OFF .. 
   | make
4. Copy the executaables to a know location on your path.
   | For example: $ mkdir ~/bin && cp -i ../src/ctmo* ~/bin/
   | Add the ~/bin to your path : 
   | in ~/.bashrc add the following line at the end:
   export PATH="$PATH:~/bin"