#!/bin/bash

if [ -z "$SLURM_JOBID" ]; then
	sbatch --nodes=32 --ntasks-per-node=1 --gid=`hostname -s` $0
else
	export OMP_NUM_THREADS=16
	srun --ntasks-per-node=1 --chdir='/gpfs/DDNgpfs3/xntmuyang/pass/run' ../bin/genlin
fi
