#!/bin/bash

cd $TRAVIS_BUILD_DIR
mkdir build_mpi
cd build_mpi
cmake -DMPI_BUILD=ON .. 
make -j 4 install
cp deps/cubature/lib/libcubature.so /usr/lib/
cp ctmo* /usr/bin/

