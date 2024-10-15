#!/bin/sh

rm -rfv build
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=/scratch1/knomura/libtorch .. && make
