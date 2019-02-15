#!/bin/bash

clang-tidy -header-filter="include/*" \
	   -p build/ src/ctmo.cpp # -fix

