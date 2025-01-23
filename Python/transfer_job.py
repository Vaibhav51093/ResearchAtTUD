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
from pyiron import ase_to_pyiron, pyiron_to_ase
from ase.io import read, write

# Project 
pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/frank/zr_vacancy_2")



q_rate = [523, 573, 623, 673, 723, 773]

#for k in range(10,12,1):
for i in q_rate:
    job = pr['nasi_2_5_random_struct_no_temp_%s'%i]
    job.transfer_from_remote()
    job.compress()
        
