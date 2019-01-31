#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    brew update 
    brew install armadillo
    brew install boost
    brew install boost-mpi
    brew install snappy 
else

    sudo update-alternatives --remove-all g++
    sudo update-alternatives --remove-all clang++
    sudo apt-get install g++-7 clang-6.0

    sudo apt-get install -y build-essential liblapack-dev cmake \
                            libarmadillo-dev libboost-all-dev libopenmpi-dev \
                            python3-numpy tree libsnappy-dev 

    sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-7 20
    sudo update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-6.0 20


fi




