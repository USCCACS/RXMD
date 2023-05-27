#!/bin/sh

CC -O2 main.cpp -L../build -L/soft/datascience/conda/2022-07-19/pytorch/torch/lib \
	-Wl,-R/soft/datascience/conda/2022-07-19/pytorch/torch/lib \
	-ltorch_cpu -ltorch_cuda -lc10 -lrxmdtorch -lstdc++
