#!/bin/sh
proj=UltrafastMat
#nodes=4356
nodes=4096
#nodes=2048
queue=default
email=knomura@usc.edu

qsub -M ${email} -A ${proj} -q ${queue} -n ${nodes} -t 1440 --mode script ./run.sh
