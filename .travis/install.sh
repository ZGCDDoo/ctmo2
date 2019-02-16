#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
    
    cd $TRAVIS_BUILD_DIR
    mkdir build
    cd build
    cmake .. && make -j 2
    cd test
    
    ./SelfConsistencyTests
    ./MarkovChainTests
    ./MarkovChainTriangle2x2Tests
    ./MarkovChainSquare2x2Tests
    ./GreenTauTests

else

    cd $TRAVIS_BUILD_DIR
    mkdir build_clang
    cd build_clang
    CXX=clang++ cmake .. && make -j 2
    make test

    cd $TRAVIS_BUILD_DIR
    mkdir build
    cd build
    cmake .. && make -j 2
    make test

fi