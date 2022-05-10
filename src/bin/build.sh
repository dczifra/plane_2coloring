#!/bin/bash

export PATH=/usr/local/cuda-10.1/bin:$PATH

time nvcc -std=c++11 -O3 -lcurand -arch=sm_35 \
    main.cu -o bin/main

echo "ECL-CUDA is ready"