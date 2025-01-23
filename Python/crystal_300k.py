# Direct lammps script implemetation in pyiron 
from pyiron import Project
import numpy as np
import pandas
from jinja2 import Template
import matplotlib.pyplot as plt 
import scipy.constants as sc
from scipy.integrate import cumtrapz
import os

# Project 
pr_1 = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/different_quench/ordered")           #! Path for postprocessing 
pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass")   

# Customizing the potentials (New potentials) [Adjust as per the Requrnments]
import pandas
padone_potential = pandas.DataFrame({
  'Name': ['NASICON_Pedone'],
  'Filename': [[]],
  'Model': ['Custom'],
  'Species': [['O', 'Na', 'Zr', 'Si', 'P']],
  'Config': [['atom_style full\n',      # Given function contains morse parametrs, i.e. it has bonds 
              '## create groups ###\n',
              'group O type 1\n',
              'group Na type 2\n',
              'group Zr type 3\n',
              'group Si type 4\n',
              'group P type 5\n',
              '\n',
              '## set charges - beside manually ###\n',
              'set group O charge -1.2000\n',
              'set group Na charge 0.6000\n',
              'set group Zr charge 2.4000\n',
              'set group Si charge 2.4000\n',
              'set group P charge 3.0000\n',
              '\n',
              'pair_style hybrid/overlay morse 15.0 mie/cut 15.0 coul/long 15.0 beck 15.0\n',
              'pair_coeff * * coul/long\n',
              'pair_coeff 1 2 beck 5.0 0 0 0 0\n',
              'pair_coeff 1 3 beck 1.0 0 0 0 0\n',
              'pair_coeff 1 4 beck 1.0 0 0 0 0\n',
              'pair_coeff 1 5 beck 1.0 0 0 0 0\n',
              'pair_coeff 1 1 beck 22.0 0 0 0 0\n',
              'pair_coeff 1 2 mie/cut 5.0 1.0 12.0 0\n',
              'pair_coeff 1 3 mie/cut 1.0 1.0 12.0 0\n',
              'pair_coeff 1 4 mie/cut 1.0 1.0 12.0 0\n',
              'pair_coeff 1 5 mie/cut 1.0 1.0 12.0 0\n',
              'pair_coeff 1 1 mie/cut 22.0 1.0 12.0 0\n',
              'pair_coeff 1 2 morse 0.023363 1.763867 3.006315\n',
              'pair_coeff 1 3 morse 0.206237 2.479675 2.436997\n',
              'pair_coeff 1 4 morse 0.340554 2.006700 2.100000\n',
              'pair_coeff 1 5 morse 0.831326 2.585833 1.800790\n',
              'pair_coeff 1 1 morse 0.042395 1.379316 3.618701\n',
              'kspace_style ewald 1.0e-8\n']]
})   


lmp_input = """\
# This script will run the calculations in one run in pyiron.
# official work of Vaibhav Arun Deshmukh, TU-Darmstadt

#------------------------------Cook and Quench method---------------------------#

#---------------------------- Atomic setup ------------------------------------#
units metal
dimension 3
boundary p p p
atom_style full
read_data structure.inp
include potential.inp
neigh_modify     delay 0
#------------------------------------------------------------------------------#

#--------------------------- Equilibrate at structure temp -----------------------------#
velocity 	    all create 600 4928459 mom yes rot yes dist gaussian
fix		        mynvt all nvt temp 300 300 0.1
timestep	    0.001
variable        dumptime  equal 100
variable        thermotime  equal 100
dump            1 all custom ${dumptime} dump.out id type xsu ysu zsu fx fy fz vx vy vz
dump_modify     1 sort id format line "%d %d %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g"
thermo_style    custom step temp pe etotal pxx pxy pxz pyy pyz pzz vol density
thermo_modify   format float %20.15g
thermo          ${thermotime}
run             60000
unfix           mynvt
#---------------------------------------------------------------------------------------------#

#-------------------------------- Equilibrate at higher temp ---------------------------------#
fix		        mynvt all nvt temp 300 300 0.1
run		        100000
unfix           mynvt
#---------------------------------------------------------------------------------------------#

#-------------------------------- Equilibrate at little lower---------------------------------#
fix		        mynvt all nvt temp 300 300 0.1
run		        100000
unfix           mynvt
#---------------------------------------- Quench NPT -----------------------------------------#
#fix		        mynpt all npt temp 300 1 0.1 iso 0.00 0.00 1.0
#fix		        mynpt all nvt temp 3000 1 0.1 
#run		        {{ step }}
#unfix           mynpt 
#--------------------------------------------------------------------------------------------#
fix		        mynpt all npt temp 300 300 0.1 iso 0.00 0.00 1.0
run		        {{ step }}
unfix           mynpt
#--------------------------------------------------------------------------------------------#
fix		        mynvt all nvt temp 300 300 0.1
run		        5000
unfix           mynvt
#--------------------------------------------------------------------------------------------#
fix		        mynve all nve 
run		        5000
unfix           mynve
#--------------------------------------------------------------------------------------------#
# Minimize the structure 
minimize 	    1.0e-4 1.0e-6 100 1000
write_data   	data.lammps
#--------------------------------------------------------------------------------------------#
"""

script_run_lmp = Template(lmp_input)

# Calculate quech rate and time 
rate = range(1, 11, 1)
#rate_1 = range(1, 11, 1)

t = [20000]     # time steps in fs 
#for i in rate:
#    time = ((4000)/i)*1000
#    t.append(round(time))

# For different colling rate 

strct_no = [1] 

for k in strct_no:
    job_mini = pr_1['boilot_623K_1_%s_ordered_mini'%k]  # Bulk Crystal
    # Get last structure    
    struct_nvt = job_mini.get_structure(iteration_step=-1)    # Last frame from minimization
    import random
    #ran = random.randrange(49284)
    for i, j in zip(t, rate):
        #ran = random.randrange(49284)
        job_name = 'boilot_443K_2_05_cryt_s_%s_1_schulze_rate_%s_k_ps'%(k,j)
        job_lmp = pr.create.job.Lammps(job_name, delete_existing_job=True)
        job_lmp.structure = struct_nvt
        job_lmp.potential = padone_potential
        job_lmp.input.control.load_string(
            script_run_lmp.render(
                step=i 
            )
        )
        job_lmp.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
        job_lmp.server.cores = 24             ### The number of cores you want for your job
        job_lmp.run()  
