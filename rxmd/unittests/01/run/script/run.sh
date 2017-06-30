#!/bin/sh

(cd init && make)
(cd src && make pre_build && make -j 24)
./rxmd
