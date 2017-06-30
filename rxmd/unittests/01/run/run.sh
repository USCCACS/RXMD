#!/bin/sh

(cd init && make)
(cd build && make pre_build && make -j 24)
./rxmd
