#!/bin/bash

nodes=4
threads=16
cores=${nodes}x${threads/16/A}

jbsub="jbsub"
# jusub+=" -interactive"
queue="x86_excl"

proj=$1
pass=$2
model=$3

mpirun="/hlt/exec/mpiwrap.sh -per-node 1"

bindir=$(readlink -f ../bin)
logdir=$(readlink -f ../log)

bin=${bindir}/${pass}

mkdir -p ${logdir}

ut=${logdir}/%J.out
err=${logdir}/%J.err

${jbsub} -queue ${queue} -proj ${proj} -name ${pass}:${model} \
	-pjobs 16 -cores ${cores} -out ${out} -err ${err} \
	OMP_NUM_THREADS=${threads} ${mpirun} ${bin} ${binopt}
