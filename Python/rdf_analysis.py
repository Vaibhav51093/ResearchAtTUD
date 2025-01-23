# The script will do the following [Using ASE and Pyiron]
# 1. RDF analysis (partial and total time avg.)
# 2. No of neighbours analysis, their indices and names 
# 3. Bond length and no of bond analysis 

# Import all the required lib
from ase import neighborlist
from ase.data import covalent_radii 
from ase.calculators.neighborlist import NeighborList
from ase.neighborlist import neighbor_list
from ase.io import read, write 
from ase import Atoms 
from scipy import sparse
import numpy as np 
from ase.io.pov import get_bondpairs 
import numpy as np
from pyiron import Project
from ase.io import read, write
from pyiron import ase_to_pyiron
import ase
import os
import time  

# look for implemetations:-  https://github.com/GardenGroupUO/GeoProps/tree/main/GeoProps
# Input file and path from pyiron 
pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/different_quench")

