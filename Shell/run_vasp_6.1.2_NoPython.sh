#!/bin/bash
module purge
module load vasp/6.1.2.NoPython
srun -n 1 vasp_std
