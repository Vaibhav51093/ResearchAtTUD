#!/bin/bash
#SBATCH -J 04d773c
#SBATCH --mail-type=END
#SBATCH --mail-user=marcus.schulze@tu-darmstadt.de
#SBATCH -e out.%x.%j
#SBATCH -o err.%x.%j
#SBATCH --mem-per-cpu=2600
#SBATCH -t 30:00:00
#SBATCH -p test7d
#SBATCH -A project01564
#SBATCH -n 24
#SBATCH --no-requeue
​
# load modules
module purge
module load intel/2021.1
module load openucx/1.10.1
module load openmpi/4.0.5
​
# print information
echo "Start Time:"
date
echo "This is Job $SLURM_JOB_ID : $SLURM_JOB_NAME on Project Number $SLURM_JOB_ACCOUNT"
echo "Nodes: $SLURM_JOB_NUM_NODES"
​
mpirun /work/home/uh84ucyh/bin/lmp_intel_cpu_openmpi -i pedone_tabular_diffusion.in -l 24core.log -sf intel
pigz -f *dump *data
​
echo "End Time:"
date
