#!/bin/bash

cppcheck -I /usr/include/armadillo_bits -I include --enable=all --inconclusive --std=posix --force ./src/ctmo.cpp

