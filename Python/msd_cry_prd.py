# Get avg. cell size and volume 
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


from random import randint
temp_list = ['473','523','573','623', '723']#,'773']                    # 623 K 
# Time = 50ps, 250ps, 37_7ps, 5ps
# NPT

job_3 = []
for temp in temp_list: 
    job_mini_3 = pr['boilot_443K_2_05_cryt_%s_NPT_diffuse'%temp]                     # crystal at 0 k after 300   
    struct_nve = job_mini_3.get_structure(iteration_step=-1)    # Last frame 
    job_name = 'boilot_443K_2_05_cryt' + "_" + temp + "_" + 'NVT_diffuse_prod'
    job_nve = pr.create_job("Lammps", job_name, delete_existing_job=True)  
    job_nve.structure = struct_nve
    job_nve.calc_md(temperature=temp, pressure=0, n_ionic_steps=5000000, n_print=100, time_step=1.0)
    job_nve.potential = padone_potential
    #job_nve.input.control['change_box'] = 'all x final '+xx_new+ ' y final ' +yy_new+' z final '+zz_new+' xy final '+xy_new+' xz final '+xz_new+' yz final '+yz_new        ### edit the box size 
    job_nve.server.queue = "normal_l2"                ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
    job_nve.server.cores = 24                         ### The number of cores you want for your job
    job_nve.run()    
    job_3.append(job_nve)