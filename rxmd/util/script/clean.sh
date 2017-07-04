#!/bin/sh

(cd init; make clean)
(cd src; make clean)
rm -fv DAT/* rxmd *.txt

find regtests/*/DAT -type f ! -name .gitignore | xargs rm -v
