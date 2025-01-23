# Direct lammps script implemetation in pyiron 
from pyiron import Project
import numpy as np
import pandas
from jinja2 import Template
import matplotlib.pyplot as plt 
import scipy.constants as sc
from scipy.integrate import cumtrapz
import os
from ase.io import read, write
from pyiron import ase_to_pyiron
from ase.io.vasp import read_vasp, write_vasp

# Project  
pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/gb/gb_111_21")   

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

#------------------------------ GB energy --------------------------------------#

#---------------------------- Atomic setup -------------------------------------#

units metal
dimension 3
boundary p p p
atom_style full
read_data structure.inp
include potential.inp
neigh_modify     delay 0

#-------------------------------------------------------------------------------#
#----------- Displace atoms and delete overlapping atoms -----------------------# 
# displace_atoms upper move 0 0 0 units lattice 
delete_atoms overlap 0.35 lower upper
 
# ---------- Define Settings ---------------------------------------------------# 
#compute csym all centro/atom fcc
compute eng all pe/atom 
compute eatoms all reduce sum c_eng 

# ---------- Run Minimization --------------------------------------------------# 
reset_timestep 0 
thermo 10 
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms 
dump 1 all cfg 25 dump.sig5_minimization_*.cfg mass type xs ys zs c_eng fx fy fz
dump_modify 1 element Al Al
min_style cg 
minimize 1e-15 1e-15 5000 5000 
undump 1

# ---------- Run Minimization 2--------------------------------------------------# 
# Now allow the box to expand/contract perpendicular to the grain boundary
reset_timestep 0 
thermo 10 
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms 
fix 1 all box/relax y 0 vmax 0.001
min_style cg 
minimize 1e-15 1e-15 5000 5000 

# ---------- Calculate GB Energy ------------------------------------------------# 
variable minimumenergy equal -14.97658984992389
variable esum equal "v_minimumenergy * count(all)" 
variable xseng equal "c_eatoms - (v_minimumenergy * count(all))" 
variable gbarea equal "lx * lz * 2" 
variable gbe equal "(c_eatoms - (v_minimumenergy * count(all)))/v_gbarea" 
variable gbemJm2 equal ${gbe}*16021.7733 
variable gbernd equal round(${gbemJm2}) 
print "GB energy is ${gbemJm2} mJ/m^2" 
 
# ---------- Dump data into Data file ------------- 
reset_timestep 0 
dump 1 all cfg 10000 dump.al_sig5_310_*.cfg mass type xs ys zs c_eng fx fy fz
dump_modify 1 element Al Al
minimize 1e-15 1e-15 ${min} ${min}
undump 1

write_restart restart.al_sig5_310_stgb

print "All done" 
#------------------------------------------------------------------------------#
"""

# Read the input script 
script_run_lmp = Template(lmp_input)
for i in range(1,2,1):
    # Import structure .cif file [Analyze the structure for its correctness seperately]
    bulk_struct = ase_to_pyiron(read(filename='111_21_nasi.cif',format='cif'))
    # Create supercel for MD simulations   
    #struct_bulk = bulk_struct.repeat([3,3,1])  
    job_name = 'boilot_gb_1_%s_ordered_mini_111_21'%i   
    job_lmp = pr_1.create.job.Lammps(job_name, delete_existing_job=True)
    #job_lmp.calc_minimize()
    job_lmp.structure = bulk_struct
    job_lmp.potential = padone_potential
    job_lmp.input.control.load_string(
        script_run_lmp.render(
            min=5000 
        )
    )
    job_lmp.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
    job_lmp.server.cores = 24             ### The number of cores you want for your job
    job_lmp.run()
    job_lmp.decompress()
