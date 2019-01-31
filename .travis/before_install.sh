#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    brew update 
    brew install armadillo
    brew install boost
    brew install boost-mpi
    brew install snappy 
else

    sudo apt-get install -y build-essential liblapack-dev cmake \
                            libarmadillo-dev libboost-all-dev libopenmpi-dev \
                            python3-numpy tree libsnappy-dev clang-6.0

fi




