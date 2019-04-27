#!/bin/bash

cd $TRAVIS_BUILD_DIR
mkdir build_mpi
cd build_mpi
CXX=mpic++ cmake -DBUILD_MPI=ON -DBUILD_TESTS=OFF ..
sudo make -j 4 install

