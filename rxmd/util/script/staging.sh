#!/bin/sh
# staging binary files with a given MD step to restart simulation

for d in DAT/??????;do
   n=`basename ${d}`
   for f in ${d}/rxff${n}-${1}.bin;do
      file=`basename ${f} | sed 's/-.*//'`
      cp -v ${f} DAT/${file}
   done
done
