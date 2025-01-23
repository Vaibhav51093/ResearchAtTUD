#!/bin/bash

#SBATCH -A p0020076
#SBATCH --array=1-11%1
#SBATCH --mail-user=deshmukh@mm.tu-darmstadt.de
#SBATCH --mail-type=REQUEUE,TIME_LIMIT
#SBATCH -e err.%j
#SBATCH -o out.%j

#SBATCH --job-name={{job_name}}
#SBATCH --chdir={{working_directory}}

{%- if run_time_max %}
#SBATCH --time={{run_time_max}}
{%- endif %}
{%- if memory_max %}
#SBATCH --mem-per-cpu={{memory_max}}
{%- endif %}
#SBATCH --ntasks={{cores}}

{{command}}
