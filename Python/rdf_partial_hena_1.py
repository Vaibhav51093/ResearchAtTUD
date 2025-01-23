from pyiron import Project, ase_to_pyiron
import matplotlib.pyplot as plt
import numpy as np
from pyiron import Project
from ase.io import read, write
from pyiron import ase_to_pyiron
import ase
import os
import time
from pyiron import Project
from pyiron import ase_to_pyiron, pyiron_to_ase
import sys  
sys.path.insert(0, '/nfshome/deshmukh/vaibhav/scripts')
import analysis_msd as ms
from scipy.optimize import curve_fit
from pyiron import Project
from pyiron import ase_to_pyiron, pyiron_to_ase
import shutil
import glob
from scipy import stats
from scipy.optimize import curve_fit
import numpy as np 
import scipy.constants as const
import pandas as pd
from scipy.optimize import curve_fit
import numpy as np 
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
import numpy as np
import os
from ase import neighborlist
from ase.data import covalent_radii 
from ase.calculators.neighborlist import NeighborList
from ase.neighborlist import neighbor_list
from ase.io import read, write 
from ase import Atoms 
from scipy import sparse
import numpy as np 
from ase.io.pov import get_bondpairs
from ovito.io import import_file
from ovito.io.ase import ovito_to_ase  

# Create an OVITO data pipeline from an external file:
#pipeline_3 = import_file("/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/glass/hena_1_glass/glass_43_623k/dump_nvt_prod.out")
pipeline_2 = import_file("/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_623k_43_hdf5/hena_1_struct_eq_big_623k_43/dump_nvt_prod.out")

# Ovito RDF analysis GB and grain 
#def setup_particle_types(frame, data):
#    types = data.particles_.particle_types_
#    types.type_by_id_(1).name = "Zr"
#    types.type_by_id_(2).name = "Si"
#    types.type_by_id_(3).name = "P"
#    types.type_by_id_(4).name = "O"
#    types.type_by_id_(5).name = "Na"
#    types.type_by_id_(4).radius = 0.3
#        
#pipeline_3.modifiers.append(setup_particle_types)
#
##pipeline_3.modifiers.append(SmoothTrajectoryModifier(window_size=10))
#
## Slice:
##pipeline_3.modifiers.append(SliceModifier(
##    distance = 0.5, 
##    normal = (0.0, 0.0, 1.0), 
##    slab_width = 10.0, 
##    inverse = True, 
##    miller = True))
#
## Insert the RDF calculation modifier into the pipeline:
#pipeline_3.modifiers.append(CoordinationAnalysisModifier(cutoff = 8.0, number_of_bins = 200, partial = True))
#
## Insert the time-averaging modifier into the pipeline, which accumulates
## the instantaneous DataTable produced by the previous modifier and computes a mean histogram.
#pipeline_3.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))
#
## Data export method 1: Convert to NumPy array and write data to a text file:
#total_rdf_crystal = pipeline_3.compute().tables['coordination-rdf[average]'].xy()
#
### Data export method 2: Use OVITO's own export function for DataTable objects:
#export_file(pipeline_3, "hena_1_big_glass_623k_partial.txt", "txt/table", key="coordination-rdf[average]")


# Ovito RDF analysis GB and grain 
def setup_particle_types(frame, data):
    types = data.particles_.particle_types_
    types.type_by_id_(1).name = "Zr"
    types.type_by_id_(2).name = "Si"
    types.type_by_id_(3).name = "P"
    types.type_by_id_(4).name = "O"
    types.type_by_id_(5).name = "Na"
    types.type_by_id_(4).radius = 0.3
        
pipeline_2.modifiers.append(setup_particle_types)

#pipeline_3.modifiers.append(SmoothTrajectoryModifier(window_size=10))

# Slice:
#pipeline_2.modifiers.append(SliceModifier(
#    distance = 0.5, 
#    normal = (0.0, 0.0, 1.0), 
#    slab_width = 10.0, 
#    miller = True))

# Insert the RDF calculation modifier into the pipeline:
pipeline_2.modifiers.append(CoordinationAnalysisModifier(cutoff = 8.0, number_of_bins = 200, partial = True))

# Insert the time-averaging modifier into the pipeline, which accumulates
# the instantaneous DataTable produced by the previous modifier and computes a mean histogram.
pipeline_2.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))

# Data export method 1: Convert to NumPy array and write data to a text file:
total_rdf_crystal = pipeline_2.compute().tables['coordination-rdf[average]'].xy()

## Data export method 2: Use OVITO's own export function for DataTable objects:
export_file(pipeline_2, "hena_1_big_cryst_623k_partial.txt", "txt/table", key="coordination-rdf[average]")
