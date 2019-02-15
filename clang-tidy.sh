#!/bin/bash

clang-tidy -header-filter="include/ctmo/MonteCarlo/*;include/ctmo/Model/*;include/ctmo/ImpuritySolver/*;include/ctom/Foundations/*" \
	   -p build/ src/ctmo.cpp # -fix

