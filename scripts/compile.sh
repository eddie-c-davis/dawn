#!/bin/bash
g++ $1.cpp -I/home/eddied/dawn/dawn/src -I/home/eddied/dawn/gtclang -I/home/eddied/dawn/gtclang/src -I/home/eddied/dawn/build/_deps/gridtools-src/include -I/home/eddied/dawn/build/_deps/googletest-src/googletest/include -std=c++11 -fopenmp -g -O3 -o $1
