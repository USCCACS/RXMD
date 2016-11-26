#!/bin/sh

(cd init; make clean)
(cd src; make clean)
rm -v DAT/*
