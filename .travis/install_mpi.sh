#!/bin/bash

cd $TRAVIS_BUILD_DIR
mkdir build_mpi
cd build_mpi
CXX=mpic++ cmake -DBUILD_MPI=ON .. 
sudo make -j 4 install

