#!/bin/bash -e
#./getStack.sh core.0 path_to_exec && cat stack.0

expectedArgs=2
if [ $# -ne $expectedArgs ]; then
  echo "Usage: $0 <corefile> <exe>"
  exit 0
fi
corefile=$1
stackfile=stack${corefile##core}
exe=$2
echo input: $corefile
echo output: $stackfile
grep -n STACK $corefile | awk -F : '{print  $1}' > lines
let s=`head -n 1 lines`+2
let f=`tail -n -1 lines`-1
sed -n ${s},${f}p $corefile | awk '{print $2}' | perl -pi -e 's/000000000/0x/g' > core.addy
addr2line -e $exe < core.addy > $stackfile
rm lines
rm core.addy
