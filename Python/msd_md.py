# MSD and diffusion coefficient calculations 
# Cook and quench method for amorphous structure ccreation
# official work by Vaibhav Deshmukh 
# Python libararies to import 
from lib2to3.pgen2.token import COLONEQUAL
from re import S
from tkinter import Y
import numpy as np
from pyiron import Project
from ase.io import read, write
from pyiron import ase_to_pyiron
import ase
import os
import time

pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/amorphous_structure") 

# Search structure:-https://icsd.fiz-karlsruhe.de/display/list.xhtml
# Use the structure after supercell low energy structure 
# Import structure .cif file [Analyze the structure for its correctness seperately]
#bulk_struct = ase_to_pyiron(read(filename='443K_nasicon_random.cif',format='cif'))                     #! Change 'file_name'                    #! Change 'file_name'

# Create supercel for MD simulations   
#struct_bulk = bulk_struct.repeat([3,3,1])  #! Choose in a such way that it should look like cube
#                                          # Also, one can build it outside and directly import the file 

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

# Cook and quench methods with diffrent cooling rates  
# 1. NVT equilibration (Timesteps=1fs, sampling and thermo(keep default))=evry 100fs, for 200K, 623K, 800K, 1400K)
from random import randint
temp_list = ['773']                    # 623 K 
# Time = 50ps, 250ps, 37_7ps, 5ps
#job_mini = pr['boilot_443K_2_05_amo_300_4_623k_1273k_new']                     # crystal at 0 k after 300   
#struct_nvt = job_mini.get_structure(iteration_step=-1)    # Last frame
#job_1 = []
#for temp in temp_list: 
#    job_name = 'boilot_443K_2_05_amo' + "_" + temp + "_" + 'NVT_diffuse'
#    job_nvt = pr.create_job("Lammps", job_name, delete_existing_job=True)  
#    job_nvt.structure = struct_nvt
#    job_nvt.calc_md(temperature=temp, pressure=None, n_ionic_steps=500000, n_print=100, time_step=1.0)
#    job_nvt.potential = padone_potential
#    job_nvt.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
#    job_nvt.server.cores = 24             ### The number of cores you want for your job
#    job_nvt.run()    
##    job_1.append(job_nvt)

#for nvt in job_1:
# Transfer all the files if the jobs are finished 
#    count = [0]
#    for x in count:
#        count.append(x+1)
#        if pr.queue_check_job_is_waiting_or_running(nvt) == True:
#            continue
#        else:
#            time.sleep(30)
#            break 
#    time.sleep(30)                             # Sleeping the process to transfer files properly 
#    job_nvt.transfer_from_remote()
#    job_nvt.compress() 
#    time.sleep(30)
#print('All the NVT calculations (amorphous) are finished............Thank you for the patience') 

# Cook and quench methods with diffrent cooling rates  
# 1. NVT equilibration (Timesteps=1fs, sampling and thermo(keep default))=evry 100fs, for 200K, 623K, 800K, 1400K)
# Time = 50ps, 250ps, 37_7ps, 5ps
# NPT
job_2 = []
for temp in temp_list:
    job_mini_2 = pr['boilot_443K_2_05_amo_%s_NVT_diffuse'%temp]                     # crystal at 0 k after 300   
    struct_npt = job_mini_2.get_structure(iteration_step=-1)    # Last frame 
    job_name = 'boilot_443K_2_05_amo' + "_" + temp + "_" + 'NPT_diffuse'
    job_npt = pr.create_job("Lammps", job_name, delete_existing_job=True)  
    job_npt.structure = struct_npt
    job_npt.calc_md(temperature=temp, pressure=1.0, n_ionic_steps=500000, n_print=100, time_step=1.0)
    job_npt.potential = padone_potential
    job_npt.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
    job_npt.server.cores = 24             ### The number of cores you want for your job
    job_npt.run()    
    job_2.append(job_npt)
# Transfer all the files if the jobs are finished 

#for npt in job_2:
#    count = [0]
#    for z in count:
#        count.append(z+1)
#        if pr.queue_check_job_is_waiting_or_running(npt) == True:
#            continue
#        else:
#            time.sleep(30)
#            break 
#    time.sleep(30)                             # Sleeping the process to transfer files properly 
#    job_npt.transfer_from_remote()
#    job_npt.compress() 
#    time.sleep(30)

#print('All the NPT calculations (amorphous) are finished............Thank you for the patience') 


#job_3 = []
#for temp in temp_list: 
#    job_mini_3 = pr['boilot_443K_2_05_amo_%s_NPT_diffuse'%temp]                     # crystal at 0 k after 300   
#    struct_nve = job_mini_3.get_structure(iteration_step=-1)    # Last frame 
#    job_name = 'boilot_443K_2_05_amo' + "_" + temp + "_" + 'NVT_diffuse_prod'
#    job_nve = pr.create_job("Lammps", job_name, delete_existing_job=True)  
#    job_nve.structure = struct_nve
#    job_nve.calc_md(temperature=temp, pressure=0, n_ionic_steps=1000000, n_print=100, time_step=1.0)
#    job_nve.potential = padone_potential
#    job_nve.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
#    job_nve.server.cores = 24             ### The number of cores you want for your job
#    job_nve.run()    
#    job_3.append(job_nve)
# Transfer all the files if the jobs are finished 

#for nve in job_3:
#    count = [0]
#    for y in count:
#        count.append(y+1)
#        if pr.queue_check_job_is_waiting_or_running(nve) == True:
#            continue
#        else:
#            time.sleep(30)
#            break 
#    time.sleep(30)                             # Sleeping the process to transfer files properly 
#    job_nve.transfer_from_remote()
#    job_nve.compress() 
#    time.sleep(30)
#    print('All the NVE calculations (amorphous) are finished............Thank you for the patience')