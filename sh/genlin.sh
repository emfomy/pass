#!/bin/bash

cores=4xA

jbsub="jbsub"
queue="x86_excl"

proj=$1
name=$2

mpirun="/hlt/exec/mpiwrap.sh -per-node 1"

bindir=$(dirname ${PWD})/bin
logdir=$(dirname ${PWD})/log

mkdir -p ${logdir}

bin=${bindir}/${name}
out=${logdir}/%J.out
err=${logdir}/%J.err

opt="-pjobs 16 -cores ${cores} -out ${out} -err ${err}"

${jbsub} -queue ${queue} -proj ${proj} -name ${name} ${opt} ${mpirun} ${bin} ${size}
