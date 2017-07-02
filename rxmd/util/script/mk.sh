#!/bin/sh

SRCDIR=`dirname $0`/../../src
cd ${SRCDIR} ; make clean; make -j 12
