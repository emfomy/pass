#!/bin/bash

cores=8xA

binopt=

jbsub="jbsub"
queue="x86_excl"

proj=$1
pass=$2
model=$3

mpirun="/hlt/exec/mpiwrap.sh -per-node 1"

bindir=$(readlink -f ../bin)
logdir=$(readlink -f ../log)

bin=${bindir}/${pass}

mkdir -p ${logdir}

out=${logdir}/%J.out
err=${logdir}/%J.err

jbopt="-pjobs 16 -cores ${cores} -out ${out} -err ${err}"

${jbsub} -queue ${queue} -proj ${proj} -name ${pass}:${model} \
	-pjobs 16 -cores ${cores} -out ${out} -err ${err} ${jbopt} \
	${mpirun} ${bin} ${binopt}
