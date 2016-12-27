#!/bin/sh

DUMPFILE=qeqdump
#DUMPFILE=rfdump

for n in 1 2; do 
   rm -fv *dump*.txt; mpirun -np ${n} ../../rxmd -o DAT-N${n}/ -in rxmd.in-N${n} && sort -n ${DUMPFILE}*.txt > log-N${n}
done
