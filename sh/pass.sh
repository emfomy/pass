#!/bin/bash

cores=8xA

jbsub="jbsub"
queue="x86_excl"

proj=$1
pass=$2
model=$3

mpirun="/hlt/exec/mpiwrap.sh -per-node 1"

bindir=$(readlink -f ../bin)
logdir=$(readlink -f ../log)

mkdir -p ${logdir}

bin=${bindir}/${pass}

for cri in AIC BIC EBIC=0.5 EBIC HDBIC HQC HDHQC ; do
	binopt="-p32 --${cri}"

	out=${logdir}/${model}_${cri}.out
	err=${logdir}/${model}_${cri}.err

	jbopt="-pjobs 16 -cores ${cores} -out ${out} -err ${err}"

	${jbsub} -queue ${queue} -proj ${proj} -name ${pass}:${model}:${cri} \
	-pjobs 16 -cores ${cores} -out ${out} -err ${err} ${jbopt} \
	${mpirun} ${bin} ${binopt}
done
