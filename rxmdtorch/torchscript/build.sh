#!/bin/sh

rm -rfv build
mkdir build
cd build
cmake ../ && make
