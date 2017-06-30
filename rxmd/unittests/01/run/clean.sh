#!/bin/sh

rm -vf rxmd DAT/* rfdump0.txt
(cd init; make clean)
(cd build; make clean)
