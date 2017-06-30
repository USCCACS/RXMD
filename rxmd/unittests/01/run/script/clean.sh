#!/bin/sh

rm -vf rxmd DAT/* rfdump0.txt
(cd init; make clean)
(cd src; make clean)
