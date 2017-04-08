#!/bin/sh

(cd init; make clean)
(cd src; make clean)
rm -fv DAT/* rxmd *.txt

find refs/*/DAT -type f | xargs rm -v
