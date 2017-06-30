#!/bin/sh

rm -vf rxmd DAT/*
(cd init; make clean)
(cd build; make clean)
