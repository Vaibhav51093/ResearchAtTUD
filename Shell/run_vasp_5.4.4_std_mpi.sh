#!/bin/bash
module purge
module load vasp/5.4.4
srun -n $1 vasp_std
