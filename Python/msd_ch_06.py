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

li = [523,573,623,673,723,773]
for i in li:
    ms.new_ovito.charge_diff(path_line='/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_2_2/minimization/hena_2_2/hena_1_struct_eq_big_%sk_6/dump_nvt_prod.out'%i,file_2='msd_na_big_ch_%sk_6.txt'%i,step=2.00, dump_write=200, step_ovito=1)