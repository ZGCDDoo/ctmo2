#!/bin/bash


cd $TRAVIS_BUILD_DIR && mkdir build && cd build && cmake .. && make -j 4
echo "Compile serial................. OK !"


cd $TRAVIS_BUILD_DIR && mkdir build_mpi && cd build_mpi && cmake -DMPI_BUILD=ON .. && make -j 4 install
echo "export PATH=\"$PATH:$HOME/bin\" "   >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH:$HOME/Installations/Libs\" "  >> ~/.bashrc
echo "Compile MPI................. OK !"


cd $TRAVIS_BUILD_DIR && cd build && make test
echo "Unit Tests................. OK !"


echo "Start of Integration Tests."
cd $TRAVIS_BUILD_DIR
travis_wait 65 python3 test/test_integration.py -v






