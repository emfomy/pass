#!/bin/bash

nodes=2
threads=16
cores=${nodes}x${threads/16/A}

# jbsub="jbsub"
jbsub="jbsub -interactive"
queue="x86_excl"

proj=$1
pass=$2
model=$3

mpirun="/hlt/exec/mpiwrap.sh -per-node 1"

bindir=$(readlink -f ../bin)
logdir=$(readlink -f ../log)

bin=${bindir}/${pass}

mkdir -p ${logdir}

for cri in HDBIC ; do
	binopt="-p1 -i128 -t10 --${cri} --prob .1 .5 .4 .5 .5"

	out=${logdir}/${model}_${cri}.out
	err=${logdir}/${model}_${cri}.err

	${jbsub} -queue ${queue} -proj ${proj} -name ${pass}:${model}:${cri} \
		-pjobs 16 -cores ${cores} -out ${out} -err ${err} \
		OMP_NUM_THREADS=${threads} ${mpirun} ${bin} ${binopt}
done
