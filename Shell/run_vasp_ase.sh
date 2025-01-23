#!/bin/bash
module purge
module use /home/groups/da_mm/modules
export MM_APPS_ROOT=/home/groups/da_mm/apps
export MM_MODULES_ROOT=/home/groups/da_mm/modules
module load vasp_vtst
export VASP_SCRIPT=/home/vd80naku/vasp/run_vasp.py
export VASP_PP_PATH=/home/vd80naku/vasp/mypps

echo $VASP_SCRIPT
echo $VASP_PP_PATH
echo $HOME
