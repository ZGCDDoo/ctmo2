#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    
    cd $TRAVIS_BUILD_DIR
    mkdir build_mpi
    cd build
    cmake -DCLANG_BUILD=ON .. && make
    make test

else

    cd $TRAVIS_BUILD_DIR
    mkdir build_clang
    cd build_clang
    cmake -DCLANG_BUILD=ON .. && make
    make test

    cd $TRAVIS_BUILD_DIR
    mkdir build
    cd build
    cmake .. && make
    make test

fi