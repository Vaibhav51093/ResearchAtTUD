# Direct lammps script implemetation in pyiron 
from tempfile import tempdir
from pyiron import Project
import numpy as np
import pandas
from jinja2 import Template
import matplotlib.pyplot as plt 
import scipy.constants as sc
from scipy.integrate import cumtrapz
import os

# Project 
pr_1 = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/amorphous_structure")                                         #! Path for postprocessing 
pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/different_quench")   

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
# This script will run the calculations in one run in pyiron-----------------------------#
# official work of Vaibhav Arun Deshmukh, TU-Darmstadt-----------------------------------#

#--------------------------- Running Diffusion Calculations -----------------------------#

#------------------------------------ Atomic setup --------------------------------------#
units metal
dimension 3
boundary p p p
atom_style full
read_data structure.inp
include potential.inp
#neigh_modify     delay 0
#-----------------------------------------------------------------------------------------#

#----------------------------------- Equilibrate NVT -------------------------------------#
reset_timestep  0
velocity 	    all create 600 49284 mom yes rot yes dist gaussian
fix		        mynvt all nvt temp {{ temp }} {{ temp }} 0.1
timestep	    0.001
variable        dumptime  equal 100
variable        thermotime  equal 100
dump            1 all custom ${dumptime} dump_nvt.out id type xsu ysu zsu fx fy fz vx vy vz
dump_modify     1 sort id format line "%d %d %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g"
thermo_style    custom step temp pe etotal pxx pxy pxz pyy pyz pzz vol
thermo_modify   format float %20.15g
thermo          ${thermotime}
run             500000
unfix           mynvt
undump          1
#-----------------------------------------------------------------------------------------#

#----------------------------------- Equilibrate NPT -------------------------------------#
reset_timestep  0
variable        dumptime  equal 100
variable        thermotime  equal 100
dump            1 all custom ${dumptime} dump_npt.out id type xsu ysu zsu fx fy fz vx vy vz
fix		        mynpt all npt temp {{ temp }} {{ temp }} 0.1 iso 1.013 1.013 1.0
thermo_style    custom step dt time atoms temp ke pe vol press lx ly lz
thermo          ${thermotime}
run             500000
unfix           mynpt
undump          1

#---------------------------------- Avg_volume/Density ------------------------------------#
##reset_timestep  0
## Getting the average volume of the system
##variable        Volume equal vol
##fix             VoluAve all ave/time 1 100000 100000 v_Volume file volume.dat
##run              100000
##reset_timestep   0

## scaling the size of the system to the average volume
##variable        sidesize equal (f_VoluAve^(1.0/3.0))    # get the volume
##variable        xlow equal xlo
##variable        ylow equal ylo
##variable        zlow equal zlo
##variable        xhig equal (xlo+${sidesize})
##variable        yhig equal (ylo+${sidesize})
##variable        zhig equal (zlo+${sidesize})
##change_box      all x final ${xlow} ${xhig} y final ${ylow} ${yhig} z final ${zlow} ${zhig}
##unfix          DensAve
##unfix           VoluAve

#------------------------------------- Production NVT -------------------------------------#
reset_timestep  0
fix		        mynvt all nvt temp {{ temp }} {{ temp }} 0.1
variable        dumptime  equal 100
variable        thermotime  equal 100
dump            1 all custom ${dumptime} dump_nvt_prod.out id type xsu ysu zsu fx fy fz vx vy vz
dump_modify     1 sort id format line "%d %d %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g"
thermo_style    custom step temp pe etotal pxx pxy pxz pyy pyz pzz vol
thermo_modify   format float %20.15g
thermo          ${thermotime}
run             1000000
unfix           mynvt
undump          1
#-------------------------------------------------------------------------------------------#
"""

script_run_lmp = Template(lmp_input)

# For different quench rate 
q_rate = [573, 773]
q_name = ['2','4','6','8']

for i in q_name:
    job_mini = pr['boilot_443K_2_05_amo_%s__k_ps'%i]   # Crystal
    struct_nvt = job_mini.get_structure(iteration_step=-1)    # Last frame from minimization
    for j in q_rate:
        job_name = 'boilot_443K_2_05_amo_' + str(i) + "_" + str(j) + '_k_ps'
        job_lmp = pr.create.job.Lammps(job_name, delete_existing_job=True)
        job_lmp.structure = struct_nvt
        job_lmp.potential = padone_potential
        job_lmp.input.control.load_string(
            script_run_lmp.render(
                temp=j 
            )
        )
        job_lmp.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
        job_lmp.server.cores = 24             ### The number of cores you want for your job
        job_lmp.run()
        job_lmp.decompress()    

# Crystal calculations 
for i in q_rate:
    job_mini = pr_1['boilot_443K_2_05_cryt_300_5_3300_300k_new']                  # Crystal minimization
    struct_nvt = job_mini.get_structure(iteration_step=-1)    # Last frame from minimization
    job_name = 'boilot_443K_2_05_cryt_%s_k_ps'%i
    job_lmp = pr.create.job.Lammps(job_name, delete_existing_job=True)
    job_lmp.structure = struct_nvt
    job_lmp.potential = padone_potential
    job_lmp.input.control.load_string(
        script_run_lmp.render(
            temp=i 
        )
    )
    job_lmp.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
    job_lmp.server.cores = 24             ### The number of cores you want for your job
    job_lmp.run()
    job_lmp.decompress()    

# Separate calculation for 10k/ps 
for i in q_rate:
    job_mini = pr_1['boilot_443K_2_05_amo_300_5_623k_1273k_new']                  # amo minimization
    struct_nvt = job_mini.get_structure(iteration_step=-1)                        # Last frame from minimization
    job_name = 'boilot_443K_2_05_amo_%s_10_k_ps'%i
    job_lmp = pr.create.job.Lammps(job_name, delete_existing_job=True)
    job_lmp.structure = struct_nvt
    job_lmp.potential = padone_potential
    job_lmp.input.control.load_string(
        script_run_lmp.render(
            temp=i 
        )
    )
    job_lmp.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
    job_lmp.server.cores = 24             ### The number of cores you want for your job
    job_lmp.run()
    job_lmp.decompress()   

print('All the Jobs are submitted and transfer them once they are finished...............!') 