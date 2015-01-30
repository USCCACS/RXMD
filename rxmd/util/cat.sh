#!/bin/sh

# default file suffix is 'pdb'. It can be specified 
# by the second argument.
if [ -z $4 ]; then
   suf=pdb
else
   suf=$4
fi

# delete all of  $suf files
rm -f *.$suf

for n0 in `seq -f %1.f $1 $2 $3`;do 
   n=`printf "%09d" ${n0}`
   echo $n -----------------------------------------

# if all files exit, copy them to working directory.
   rm -f temp 
   for dir in `seq -f %06g 0 59`;do
       cat DAT/${dir}/rxff${dir}-${n}.${suf} >> temp 
   done

   if [ $suf == "pdb" ]; then
      sort -n -k 4 temp > $n.$suf 
      echo END >> $n.$suf 
   fi

done

# concatinate them into one file
cat *.$suf > all
mv all all.pdb
rm -f temp
