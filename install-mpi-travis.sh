#!/bin/bash

cd $TRAVIS_BUILD_DIR
mkdir build_mpi
cd build_mpi
cmake -DMPI_BUILD=ON .. 
make -j 4 install
echo "export PATH=\"$PATH:$HOME/bin\" "   >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH:$HOME/Installations/Libs\" "  >> ~/.bashrc

