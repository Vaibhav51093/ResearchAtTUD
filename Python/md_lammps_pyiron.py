#-------------------------------------------MD_LAMMPS_RUN_SCRIPT_USING_PYIRON--------------------------------------------#
#----------------------------------------------------------------------------Created_By_Vaibhav_Arun_Deshmukh------------#
#---Instructions---------------------------------------------------------------------------------------------------------#
# Ensure the following packages before running your script:-                                                             #
# 1. Pyiron, ASE package and lammps configuration with                                                                   #
#    local PC and HPC.                                                                                                   #
# 2. Check the needed lammps subpackages.                                                                                #     
#                                                                                                                        #
# Note:- 1. #! comment indicates lammps input parameters                                                                 #
#           and they can be altered as per the requirnments.                                                             #
#        2. Create the unitcell seperately and just import                                                               #
#           .cif file and prepare supercell inside.                                                                      #
#        3. This is just an processing script and later the                                                              #
#           results can be analysze using post processing                                                                #
#           script. Therefore it is important to remember                                                                #
#           project addresse to import the results.                                                                      #
#------------------------------------------------------------------------------------------------------------------------#
# Python libararies to import 
from ast import Break
from this import d
import numpy as np
from pyiron import Project
from ase.io import read, write
from pyiron import ase_to_pyiron
import ase
import os
import time

# Define the parameters for the simulations (Mini, NVT, NPT) [I will write it letter]
# Purpose:- To define intially all the things and then you dont need to look at the 
#           detailed script  
# A. Minimization
# B. NVT
# C. NPT 


pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/crystalline")           #! Path for postprocessing 

# Import structure .cif file [Analyze the structure for its correctness seperately]
bulk_struct = ase_to_pyiron(read(filename='443K_nasicon_random.cif',format='cif'))                     #! Change 'file_name'

# Create supercel for MD simulations   
struct_bulk = bulk_struct.repeat([3,3,1])  #! Choose in a such way that it should look like cube
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

# 1. Energy Minimization 
# Create job name 
job_name = 'boilot_443K_2_05_mini'
job_mini = pr.create.job.Lammps(job_name, delete_existing_job=True)

# Specifing the bulk structure to the job 
job_mini.structure = struct_bulk

# Modifying the input strcuture 
job_mini.calc_minimize(e_tol=1.0e-4, f_tol=0.0001, max_iter=1.0e-6, n_print=10, style='cg', pressure=None)
job_mini.potential = padone_potential
#job_mini.input.control['min_style'] = 'cg' # This is how we can edit lammps input file 
#print(job_mini.input.control)
job_mini.server.queue = "normal_l2"        ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
job_mini.server.cores = 24                 ### The number of cores you want for your job

# Run minimization
job_mini.run()
# I have defined the following as new class in generic.py and queuestatus.py files in pyiron
# path for generic.py:- /nfshome/deshmukh/miniconda3/envs/vaibhav/lib/python3.10/site-packages/pyiron_base/project/generic.py
# path for queuestatus.py:- /nfshome/deshmukh/miniconda3/envs/vaibhav/lib/python3.10/site-packages/pyiron_base/server/queuestatus.py  
count = [0]
for x in count:
  count.append(x+1)
  if pr.queue_check_job_is_waiting_or_running(job_mini) == True:
    continue
  else:
    break 
time.sleep(10)                                            # Sleeping the process to transfer files properly 
#job_mini.status.collect=True
job_mini.transfer_from_remote()
job_mini.compress() 
print('Transfered' + ' ' + job_name)
time.sleep(5)
struct_nvt = job_mini.get_structure(iteration_step=-1)    # Last frame from minimization

# 2. NVT equilibration (Timesteps=1fs, sampling and thermo(keep default))=evry 100fs, for 200K, 623K, 800K, 1400K)
temp_list = ['200','623','800','1400']
#job_nvt_id = []
job = []
for temp in temp_list: 
  job_name = 'boilot_443K_2_05_NVT' + "_" + temp
  job_nvt = pr.create_job("Lammps", job_name, delete_existing_job=True)
  job_nvt.structure = struct_nvt
  #job_nvt_id.append(job_nvt)
  job_nvt.calc_md(temperature=temp, pressure=None, n_ionic_steps=1000, n_print=100, time_step=1.0)
  job_nvt.potential = padone_potential
  #print(job_mini.input.control)
  job_nvt.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
  job_nvt.server.cores = 24             ### The number of cores you want for your job
  job_nvt.run()    
  #job_nvt_id.append(job_nvt.queue_id)
  job.append(job_nvt)
  # Check input file for NVT 
  #print(job_nvt.input.control)

# Transfer all the files if the jobs are finished 
for nvt in job:
  count = [0]
  for x in count:
    count.append(x+1)
    if pr.queue_check_job_is_waiting_or_running(nvt) == True:
      continue
    else:
      break 
  time.sleep(10)                             # Sleeping the process to transfer files properly 
  #job_nvt.status.collect=True
  nvt.transfer_from_remote()
  nvt.compress() 
  print('Transfered NVT calculations')
  
# 3. NPT equilibration 
press_list = ['1.00']     # list of different pressure  
job_1 = []
for temp in temp_list: 
  job_out_nvt = pr['boilot_443K_2_05_NVT' + "_" + temp]
  job_name = 'boilot_443K_2_05_NPT' + "_" + temp
  job_npt = pr.create_job("Lammps", job_name, delete_existing_job=True)
  job_npt.structure = job_out_nvt.get_structure(iteration_step=-1)
  job_npt.calc_md(temperature=temp, pressure=press_list, n_ionic_steps=1000, n_print=100, time_step=1.0)
  job_npt.potential = padone_potential
  #print(job_mini.input.control)
  job_npt.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
  job_npt.server.cores = 24             ### The number of cores you want for your job
  job_npt.run()    
  # Check input file for NPT 
  #print(job_npt.input.control)
  job_1.append(job_npt)

#Transfer all files once the jobs are finished 
for npt in job_1:
  count = [0]
  for x in count:
    count.append(x+1)
    if pr.queue_check_job_is_waiting_or_running(npt) == True:
      continue
    else:
      break 
  time.sleep(10)                             # Sleeping the process to transfer files properly 
  #job_nvt.status.collect=True
  npt.transfer_from_remote()
  npt.compress() 
  print('Transfered NPT calculations')

print('The Calculations are finished')
