#!/bin/bash
#SBATCH -J pure_Ni_sf_path
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deshmukh@mm.tu-darmstadt.de
#SBATCH -e ./err_files/err.%x.%j
#SBATCH -o ./out_files/out.%x.%j
#SBATCH -t 2:00:00
#SBATCH -p devel
#SBATCH -A lowsfehea
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
​
# load modules
module purge
#module load Stages/2020
#module load Intel/2021.2.0-GCC-10.3.0
#module load OpenMPI/4.1.1
#module load mpi4py/3.0.3-Python-3.8.5
#module load GCCcore/.10.3.0
#module load Python/3.8.5
#module load SciPy-Stack/2021-Python-3.8.5
​
module load Intel/2021.4.0
module load OpenMPI/4.1.1
module load mpi4py/3.1.3
module load Python/3.9.6
module load SciPy-Stack/2021b
​
# Setting/updating environmental variables
export PYTHONPATH="${PYTHONPATH}:/p/home/jusers/deshmukh/juwels/Some_useful_codes/useful_common_classes/Atomistic_configurations/3D_systems_of_different_geometries/"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/p/home/jusers/nag1/juwels/Softwares/LAMMPS_installations/intel_2021_2_0_GCC_10_3_0_openmpi_4_1_1/lammps-29Oct20/build"
export PYTHONPATH="${PYTHONPATH}:/p/home/jusers/nag1/juwels/Softwares/LAMMPS_installations/intel_2021_2_0_GCC_10_3_0_openmpi_4_1_1/lammps-29Oct20/python"
export lammps_exe="/p/home/jusers/nag1/juwels/Softwares/LAMMPS_installations/intel_2021_2_0_GCC_10_3_0_openmpi_4_1_1/lammps-29Oct20/build/lmp_intel_cpu_openmpi"
​
​
# print information
echo "Start Time:"
date
echo "This is Job $SLURM_JOB_ID : $SLURM_JOB_NAME on Project Number $SLURM_JOB_ACCOUNT"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "PYTHONPATH: $PYTHONPATH"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
​
srun python3 $HOME/Some_useful_codes/Specific_problems/Path_to_stacking_fault_creation/sf_creation_path_in_fcc_using_fix_atoms_in_plane_constraint.py ./pure_Ni_sf_path.config
​
echo "End Time:"
date