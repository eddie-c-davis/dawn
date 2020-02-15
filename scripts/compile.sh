#!/bin/bash
g++ $1.cpp -I/home/eddied/Work/dawn/dawn/src -I/home/eddied/Work/dawn/gtclang -I/home/eddied/Work/dawn/build/_deps/gridtools-src/include -I/home/eddied/Work/dawn/build/_deps/googletest-src/googletest/include -c -o $1.o
