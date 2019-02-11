#!/bin/bash

clang-tidy -header-filter=./src/* -p build/ src/ctmo.cpp -fix

