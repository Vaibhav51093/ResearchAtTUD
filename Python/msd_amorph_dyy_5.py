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

li =  [523,573,623,673,723,773]
#li = [623]
no = [5,5,5,5,5,5]

# Run over several measurements

input_files = ["/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass/boilot_623K_2_05_amorph_523_k_ps_new_5_hdf5/boilot_623K_2_05_amorph_523_k_ps_new_5/dump_nvt_prod.out",
               "/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass/boilot_623K_2_05_amorph_573_k_ps_new_5_hdf5/boilot_623K_2_05_amorph_573_k_ps_new_5/dump_nvt_prod.out",
               "/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass/boilot_623K_2_05_amorph_623_k_ps_new_5_hdf5/boilot_623K_2_05_amorph_623_k_ps_new_5/dump_nvt_prod.out",
               "/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass/boilot_623K_2_05_amorph_673_k_ps_new_5_hdf5/boilot_623K_2_05_amorph_673_k_ps_new_5/dump_nvt_prod.out",
               "/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass/boilot_623K_2_05_amorph_723_k_ps_new_5_hdf5/boilot_623K_2_05_amorph_723_k_ps_new_5/dump_nvt_prod.out",
               "/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass/boilot_623K_2_05_amorph_773_k_ps_new_5_hdf5/boilot_623K_2_05_amorph_773_k_ps_new_5/dump_nvt_prod.out",]

for i,j,k in zip(li,input_files,no):
    print(j)
    m.ovito_msd.diff_3(path_line=j,file_2='msd_na_amo_dyy_%s_%s.txt'%(i,k),step=1.0, dump_write=100, step_ovito=10)
