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
import warnings
warnings.filterwarnings("ignore")

pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/von_alpen")

pr.update_from_remote()

print('Done-------')
