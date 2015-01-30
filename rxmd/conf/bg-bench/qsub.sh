#!/bin/sh
nodes=1
mode=c1
proj=PetaSimNano
queue=default
output=rxmd-${nodes}-${mode}

if [ -z $1 ]; then
  dep=""
else
  dep="--dependencies $1"
fi

qsub -t 30 ${dep} -q ${queue} -n ${nodes} -A ${proj} --mode ${mode} -O ${output} ./rxmd
