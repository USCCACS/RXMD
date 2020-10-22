#!/bin/sh

for step in `seq 0 3`;do 
for f in 1_1_1 2_2_1 2_2_4 4_4_4 8_8_4 8_8_16 16_16_16 32_32_16 32_32_64 64_64_64 66_66_64; do 
  vp=`echo $f | sed -e s'/.*-//' -e 's/_/ /g'`
  np=`echo $f | sed -e s'/.*-//' -e 's/_/\*/g' | bc`
  echo "=== $f    `date` ==="
  aprun -n ${np} ./rxmd --vprocs ${vp} --short_rep shortrep.in --random_velocity | tee log-${step}-${np}
done
done

