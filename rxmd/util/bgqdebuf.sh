#!/bin/sh

for f in core.*;do  echo -n "$f  " ; 
  n=`grep 'While executing instruction at' ${f}  | sed 's/.*\.//' `
  addr2line -e rxmd ${n}
done
