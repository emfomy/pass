#!/bin/bash

nodes=8
threads=16
cores=${nodes}x${threads/16/A}

jbsub="jbsub"
queue="x86_excl"

proj=$1
pass=$2
model=$3

binopt="-p32 --HDBIC --prob .1 .5 .4 .5 .5"

mpirun="/hlt/exec/mpiwrap.sh -per-node 1"

bindir=$(readlink -f ../bin)
logdir=$(readlink -f ../log)

bin=${bindir}/${pass}

mkdir -p ${logdir}

out=${logdir}/${model}.out
err=${logdir}/${model}.err

${jbsub} -queue ${queue} -proj ${proj} -name ${pass}:${model} \
	-pjobs 16 -cores ${cores} -out ${out} -err ${err} \
	OMP_NUM_THREADS=${threads} ${mpirun} ${bin} ${binopt}
