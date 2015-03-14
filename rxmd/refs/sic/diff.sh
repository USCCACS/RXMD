#!/bin/sh

echo "=== DAT directory difference ===: $1 $2"
diff -rq $1/DAT $2/DAT

echo "=== stdout difference ===: $1 $2"
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $1/std.out  > tmp1
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $2/std.out  > tmp2
diff tmp1 tmp2
rm tmp1 tmp2
