#!/bin/sh

icpx -std=c++17 -g -fiopenmp -fopenmp-targets=spir64 -D__STRICT_ANSI__ nn.cpp && iprof ./a.out nh3.xyz Pb.20t-20t.nn-00005_asc
