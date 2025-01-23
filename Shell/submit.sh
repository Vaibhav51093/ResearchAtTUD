#!/bin/sh
#SBATCH -A project01973 #special00003 
#SBATCH --job-name=diagonal_less 
#SBATCH --mail-user=kim@mm.tu-darmstadt.de 
#SBATCH --mail-type=END
#SBATCH --nodes=1 
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=96 
#SBATCH --cpus-per-task=1
#SBATCH --output=out.%j
#SBATCH --error=err.%j
#SBATCH --time=24:00:00 
#SBATCH --exclusive
#SBATCH --mem-per-cpu=3500
##SBATCH -n 2

##module load lammps/2021.10.27_ml
#module purge
#module load INTEL/toolchain

module purge
ml intel/2021.1 openmpi/4.1.1


# print information
echo "Start Time:"
date
echo "This is Job $SLURM_JOB_ID : $SLURM_JOB_NAME on Project Number $SLURM_JOB_ACCOUNT"
echo "Nodes: $SLURM_JOB_NUM_NODES"

srun /home/vd80naku/lammps/lammps/src/lmp_intel_cpu_openmpi -i Diffusion.in -l logfile.log -sf intel

echo "End Time:"
date
