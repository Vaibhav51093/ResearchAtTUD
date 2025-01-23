#!/bin/bash
module purge
module load mkl/2021.2.0
module load mpi/2021.2.0
exec vasp_std
module purge
