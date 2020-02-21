#!/bin/bash
DAWN_PATH=/home/eddied/dawn
nvcc $1.cu -DGRIDTOOLS_DAWN_CUDA -DOPTBACKEND=cuda -arch=sm_60 -g -O3 -std=c++11 -I${DAWN_PATH}/dawn/src -I/${DAWN_PATH}/gtclang -I${DAWN_PATH}/gtclang/src -I${DAWN_PATH}/build/_deps/gridtools-src/include -o ${1}_cu
