#!/bin/bash
#SBATCH -A p0020075
#SBATCH --array=1-11%1
#SBATCH --job-name=AIMD
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --mail-user= XXXXX
#SBATCH --mail-type= XXXXX
#SBATCH --cpus-per-task=1
#SBATCH --output=out.%A_%a
#SBATCH --error=err.%A_%a
#SBATCH --mem-per-cpu=3500
#SBATCH --time=24:00:00
#SBATCH --no-requeue


module purge
module load intel/2022.1
module load intelmpi/2022.1
module load fftw/3.3.8

current_run=$(printf "%02d" $SLURM_ARRAY_TASK_ID)    
n=$(($SLURM_ARRAY_TASK_ID-1))                        
previous_run=$(printf "%02d" $n)                     


# *** start of job script ***
# Note: The current working directory at this point is
# the directory where sbatch was executed.



### CONTCAR to POSCAR + backup of files; except of the very first run
if (( $SLURM_ARRAY_TASK_ID > $SLURM_ARRAY_TASK_MIN ))
then
    cp CONTCAR  POSCAR
    cp CONTCAR  CONTCAR_$previous_run
    cp OUTCAR   OUTCAR_$previous_run
    cp XDATCAR  XDATCAR_$previous_run
    rm CONTCAR XDATCAR OUTCAR
fi


### Backup of POSCAR + start of calculation; except of the last tidy-up run     
if (( $SLURM_ARRAY_TASK_ID < $SLURM_ARRAY_TASK_MAX ))                           
then                                                                            
    cp POSCAR  POSCAR_$current_run
    srun /home/groups/da_mm/codes/vasp/vasp.6.3.0+VTST/bin/vasp_gam
fi


### Tidy up in last run
if (( $SLURM_ARRAY_TASK_ID == $SLURM_ARRAY_TASK_MAX ))
then
    rm XDATCAR OUTCAR POSCAR CHG CHGCAR vasprun.xml PROCAR DOSCAR REPORT EIGENVAL WAVECAR
    gzip OUTCAR_*
fi