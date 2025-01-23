#!/bin/bash
module purge
#module load mkl/2021.2.0
#module load mpi/2021.2.0

# Copy kernal file for vasp VDW calculations 
cp /home/vd80naku/pyiron/resources/vasp/potentials/vdw_kernel.bindat ./ 

#module load vasp/6.1.2
module load intel/2022.1
module load intelmpi/2022.1
module load fftw/3.3.10

module load vasp_vtst/6.3.0 

srun -n $1 vasp_std
