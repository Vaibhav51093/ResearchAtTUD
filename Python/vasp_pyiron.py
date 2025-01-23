# Script that will perform K space calculation convergence 
from ast import Break
from this import d
import numpy as np
from pyiron import Project
from ase.io import read, write
from pyiron import ase_to_pyiron
import ase
import os
import time
  
pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/structure/lammps_vasp")

# 1. Kspace analysis 
job = []
k_spacings = [0.15, 0.13, 0.11, 0.09, 0.07]
job_mini = pr['boilot_443K_2_05_mini_1']                    # Loading each job 
struct_mini = job_mini.get_structure(iteration_step=-1)     # Last frame from minimization
for d in k_spacings:
    j = pr.create.job.Vasp(f"KMesh_{d}".replace(".", "_"))
    j.structure = struct_mini
    j.set_kpoints(k_mesh_spacing=d, scheme="GC")
    j.set_occupancy_smearing(smearing="gaussian")
    j.input.incar["ISIF"] = 3
    j.input.incar["KPAR"] = 6
    j.input.incar["NCORE"] = 4
    j.input.incar["EDIFF"] = 1e-8
    j.input.incar["ALGO"] = "Normal"
    j.input.incar["ENCAT"] = 600
    j.server.cores = 24
    j.server.queue = "normal_l2"
    j.run()
    job.append(j)

for nvt in job:
  count = [0]
  for x in count:
    count.append(x+1)
    if pr.queue_check_job_is_waiting_or_running(nvt) == True:
        continue
    else:
        time.sleep(30)
        break 
  time.sleep(30)                             # Sleeping the process to transfer files properly 
  #job_nvt.status.collect=True
  nvt.transfer_from_remote()
  nvt.compress() 
  time.sleep(5)
  print('Transfered structure')
print('All vasp structured transfered')