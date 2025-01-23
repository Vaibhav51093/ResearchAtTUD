#!/bin/bash

#SBATCH -A project01973
#SBATCH --mail-user=deshmukh@mm.tu-darmstadt.de
#SBATCH --mail-type=REQUEUE,TIME_LIMIT
#SBATCH -e err.%j
#SBATCH -o out.%j

#SBATCH --job-name=pi_1
#SBATCH --chdir=/work/scratch/vd80naku/pyiron/projects/remote/NASICON/project/padone_pot/glass_latest/replica/nasi_1_5_random/523
#SBATCH --time=1440
#SBATCH --mem-per-cpu=3800
#SBATCH --ntasks=96

# load modules
module purge
ml intel/2021.1 openmpi/4.1.1

# Template LAMMPS script
template_script="control.inp"

# Run LAMMPS script five times
for run in {1..5}
do
    echo "Running LAMMPS - Run ${run}"
    echo "Vaibhav's official work."	
    # Generate a copy of the template script with a unique name
    script_name="lammps_script_${run}.in"
    cp "${template_script}" "${script_name}"

    # Replace the file names in the copied script
    sed -i "s/log_n.lammps/log_n_${run}.lammps/" "${script_name}"
    sed -i "s/minimized.data/minimized_${run}.data/" "${script_name}"
    sed -i "s/restart.data/restart_${run}.data/" "${script_name}"
    sed -i "s/data_n.lammps/data_n_${run}.lammps/" "${script_name}"
    sed -i "s/dump_nvt_prod.out/dump_nvt_prod_${run}.out/" "${script_name}"

    # Run LAMMPS with random seed from Bash script
    srun /home/vd80naku/lammps/lammps/src/lmp_intel_cpu_openmpi -i "${script_name}" -v random $RANDOM 
    # Execute the LAMMPS script
    #lmp -i "${script_name}"

    echo "LAMMPS Run ${run} completed"
done

