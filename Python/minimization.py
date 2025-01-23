# Minimizing all the relavent structures to figure out the effect 
from lib2to3.pgen2.token import COLONEQUAL
import numpy as np
from pyiron import Project
from ase.io import read, write
from pyiron import ase_to_pyiron
import ase
import os
import time

pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/amorphous_structure")           #! Path for postprocessing 

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


job = []
list = ['1','2','3']
iter = 0
#job_mini_1 = pr['boilot_443K_2_05_cry6t_23_4']  # Crystal
job_mini_2 = pr['boilot_443K_2_05_amo_300_4_623k_1273k_new']
# Minimizing all the structures 
for i in list:
    #if iter == 0:
    #    job_name = 'boilot_443K_2_05_cryt_623_5'
    #    job_mini = pr.create.job.Lammps(job_name, delete_existing_job=True)
    #    job_mini.structure = job_mini_1.get_structure(iteration_step=-1)
    #    job_mini.calc_minimize()
    #    job_mini.potential = padone_potential
    #    job.append(job_name)
    if iter == 0:
        job_name = 'boilot_443K_2_05_amo_300_5_623k_1273k_new'
        job_mini = pr.create.job.Lammps(job_name, delete_existing_job=True)
        job_mini.structure = job_mini_2.get_structure(iteration_step=-1)
        job_mini.calc_minimize()
        job_mini.potential = padone_potential
        job.append(job_mini)
    else:
        break
    iter += 1
    job_mini.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
    job_mini.server.cores = 24             ### The number of cores you want for your job
    job_mini.run()    

    count = [0]
    for x in count:
        count.append(x+1)
        if pr.queue_check_job_is_waiting_or_running(job_mini) == True:
            continue
        else:
            time.sleep(30)
            break 
    time.sleep(30)                                # Sleeping the process to transfer files properly 
    job_mini.transfer_from_remote()
    job_mini.compress() 
    time.sleep(5)
    print('Transfered NVT calculations')