#!/bin/sh

(cd init; make)
(cd build; make)
./rxmd
