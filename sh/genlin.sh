#!/bin/bash

if [ -z "$SLURM_JOBID" ]; then
	sbatch --nodes=1 --ntasks-per-node=1 --qos=umax-32 --time=60 --gid=`hostname -s` $0
else
	export OMP_NUM_THREADS=32
	srun --ntasks-per-node=1 --chdir='$(pwd)' ../bin/genlin
fi
