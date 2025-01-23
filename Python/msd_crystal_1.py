#!/usr/bin/env python
# coding: utf-8

# # Form for Calculation of Diffusivity

# ## Reading of dump-files to generate txt.file containing MSD (tracer diffusivity), netMSD (ionic diffusivity), MSD-com (center of mass corrected), netMSD-com: (each for x,y,z)  
# ### ca. 30 minutes

from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
import numpy as np
import os
#from pyiron import Project, ase_to_pyiron
import matplotlib.pyplot as plt
import numpy as np
#from pyiron import Project
from ase.io import read, write
#from pyiron import ase_to_pyiron
import ase
import os
import time
import sys

sys.path.append(os.path.abspath("/nfshome/deshmukh/vaibhav/scripts"))
import analysis_msd as m 

li =  [523] #,623,673,723,773]

# Run over several measurements

input_files = ["/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass_1/boilot_623K_2_05_cryt_523_k_ps_new_3_hdf5/boilot_623K_2_05_cryt_523_k_ps_new_3/dump_nvt_prod.out"]

for i,j in zip(li,input_files):
    print(j)
    m.ovito_msd.diff_2(path_line=j,file_2='msd_na_%s_3.txt'%i,step=1.0, dump_write=100, step_ovito=10)
