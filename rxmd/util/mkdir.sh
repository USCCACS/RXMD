#!/bin/sh

if ! [ -d ./DAT ]; then
   echo "run this script from working directory."
   exit
fi

for d in DAT/rxff*; do 
  dir=`echo $d | sed 's/.*rxff//'`
  mkdir -v DAT/${dir}
done
