#!/bin/bash

cppcheck -I /usr/include/armadillo_bits --enable=all --inconclusive --std=posix --force ./src/ctmo.cpp 

