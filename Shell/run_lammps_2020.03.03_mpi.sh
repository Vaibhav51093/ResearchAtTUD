#!/bin/bash

# load modules
module purge
ml intel/2021.1 openmpi/4.1.1 

# run the simulation
srun /home/vd80naku/lammps/lammps/src/lmp_intel_cpu_openmpi -i control.inp

