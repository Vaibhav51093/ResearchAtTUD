#!/bin/bash

#SBATCH -A p0020646  #p0020076 #project01973  #p0020076 #(Marcel) #project01973 #(mine) #project01964 #(sun)
#SBATCH --mail-user=deshmukh@mm.tu-darmstadt.de
#SBATCH --mail-type=All
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
