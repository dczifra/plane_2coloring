#!/bin/bash

export PATH=/usr/local/cuda-10.1/bin:$PATH

time nvcc -std=c++11 -g -lcurand -arch=sm_35 \
    cpp/main.cu -o bin/main

echo "ECL-CUDA is ready"