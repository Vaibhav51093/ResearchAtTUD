#!/bin/bash

#SBATCH -A p0020646 #project01973  #p0020076  #project01973  #p0020646  #p0020076 #project01973  #p0020076 #(Marcel) #project01973 #(mine) #project01964 #(sun)
#SBATCH --mail-user=deshmukh@mm.tu-darmstadt.de
#SBATCH --mail-type=All
#SBATCH -e err.%j
#SBATCH -o out.%j

#SBATCH --job-name=pi_32515
#SBATCH --chdir=/work/scratch/vd80naku/pyiron/projects/remote/NASICON/new_glass/na2sio3_amorphous_large_lmp_hdf5/na2sio3_amorphous_large_lmp
#SBATCH --time=30
#SBATCH --mem-per-cpu=3800
#SBATCH --ntasks=96

python -m pyiron_base.cli wrapper -p /work/scratch/vd80naku/pyiron/projects/remote/NASICON/new_glass/na2sio3_amorphous_large_lmp_hdf5/na2sio3_amorphous_large_lmp -f /work/scratch/vd80naku/pyiron/projects/remote/NASICON/new_glass/na2sio3_amorphous_large_lmp.h5/na2sio3_amorphous_large_lmp