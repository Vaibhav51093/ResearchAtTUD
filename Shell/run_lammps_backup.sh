#!/bin/bash
#module load lammps/2020.03.03

#srun lmp -in control.inp
#mpiexec -n $1 --oversubscribe lmp_mpi -in control.inp;

# load modules
#module purge
#ml intel/2020.4 openmpi/4.1.1 

# run the simulation
#srun --mpi=openmpi /home/vd80naku/lammps/lammps/src/Obj_intel_cpu_openmpi -i control.inp
# load modules
module purge
#ml intel/2020.4 openmpi/4.1.1
#ml intel/2021.1 openucx/1.11.2 intelmpi/2021.1
#ml intel/2021.1 openucx/1.12.0 intelmpi/2021.1
#ml intel/2021.1 intelmpi/2021.1
ml intel/2021.1 openmpi/4.1.1 


# run the simulation
#srun --mpi=openmpi /home/vd80naku/lammps/lammps/src/Obj_intel_cpu_openmpi -i control.inp
srun /home/vd80naku/lammps/lammps/src/lmp_intel_cpu_openmpi -i control.inp
#srun --mpi=intelmpi /home/vd80naku/lammps/lammps/src/lmp_intel_cpu_intelmpi -i control.inp

