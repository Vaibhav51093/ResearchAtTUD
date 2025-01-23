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
import sys
sys.path.insert(0, '/nfshome/deshmukh/vaibhav/scripts')
import analysis_msd as ms

# Project 
pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization")

q_rate = [43]#, 47, 4, 13, 42] #, 10, 5, 26, 28, 30, 25] #42
q_tem = [523, 573, 623, 673, 723, 773]  

for i in q_rate:
    for j in q_tem: 
        ms.ovito_msd.diff_ch(job=pr['hena_1_struct_eq_big_%sk_%s'%(j,i)],filename='dump_nvt_prod.out',file_2='msd_na_big_ch_%sk_%s.txt'%(j,i),step=2,dump_write=200,step_ovito=10)
