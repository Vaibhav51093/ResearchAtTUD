#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Libraries I need to calculate  1. MSD, 2. Diffusion coefficient and 3. Ionic conductivity
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pylab as plt
from pyiron.project import Project
import ase.units as units
import pandas
from asap3 import*     # to plot RDF remember the function 
from asap3.analysis.rdf import RadialDistributionFunction 
from pyiron import pyiron_to_ase
from ase.io import read, write
from ase.io.trajectory import Trajectory
import matplotlib as mpl 
from pylab import cm
from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii


# In[23]:


# Creating the trajectory first and later we will import in polypy for analysis 
# Import projects and load the structure to analyze 
pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/amorphous_structure")     # Project name
job_loaded_1 = pr['boilot_443K_2_05_amo_300_4_623k_1273k_new']
job_loaded_1.decompress()
os.listdir(job_loaded_1.working_directory)


# In[37]:


dump_file = job_loaded_1['dump.out']
print(dump_file)


# In[34]:


from ase.io.lammpsrun import read_lammps_dump_text, read_lammps_dump_binary 
data = read_lammps_dump_binary(dump_file)


# In[16]:


from ase.io.trajectory import Trajectory
xyz = []
for i in range(1,51,1):
    struct = pyiron_to_ase(job_loaded_1.get_structure(iteration_step=i))
    xyz.append(struct)
write('3300_300k_new.xyz', xyz)
a = read('3300_300k_new.xyz')
write('HISTORY', a, format='dlp4')


# In[17]:


from polypy import read_dl_poly as dlpy
from polypy import read as rd
from polypy.msd import msd as msd
#from polypy.msd import RegionalMSD
from polypy import analysis
from polypy import utils as ut
from polypy import plotting


# In[18]:


data = rd.read_history("HISTORY", ["Na"])


# In[21]:


positions = job_loaded_1['output/generic/positions']
steps = job_loaded_1['output/generic/steps']
steps


# In[22]:


trajectory = ase.io.read('3300_300k_new.pdb')


# In[23]:


import mdtraj as md
md.load('3300_300k_new.pdb')


# In[27]:


import os
import json
from pymatgen.core import Structure
from pymatgen.analysis.diffusion_analyzer import DiffusionAnalyzer, get_arrhenius_plot, get_extrapolated_conductivity
from pymatgen_diffusion.aimd.pathway import ProbabilityDensityAnalysis


# In[33]:


traject = Structure.from_file("3300_300k_new.xyz")

