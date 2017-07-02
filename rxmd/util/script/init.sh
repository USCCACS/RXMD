#!/bin/sh

INITDIR=`dirname $0`/../../init
cd ${INITDIR}; make clean; make -j 12
