import argparse
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
import numpy as np
from scipy import integrate
from ovito.io.ase import ase_to_ovito, ovito_to_ase
from scipy.constants import codata
from ase.io import write
from ase.io.lammpsdata import write_lammps_data

# Setup command line argument parsing
parser = argparse.ArgumentParser(description="Convert LAMMPS dump to data file.")
parser.add_argument("filename", help="Path to the LAMMPS dump file.")
args = parser.parse_args()

# Use the filename from the command line
filename = args.filename

def setup_particle_types(frame, data):
    types = data.particles_.particle_types_
    types.type_by_id_(1).name = "O"
    types.type_by_id_(2).name = "Na"
    types.type_by_id_(3).name = "Zr"
    types.type_by_id_(4).name = "Si"
    types.type_by_id_(5).name = "P"
    #types.type_by_id_(5).radius = 0.3

# Load the dump file
pipeline = import_file(filename, multiple_frames=True)

# Add modifiers
pipeline.modifiers.append(setup_particle_types)
pipeline.modifiers.append(WrapPeriodicImagesModifier())

# Compute the pipeline and convert to ASE atoms object
data = pipeline.compute()
atom = ovito_to_ase(data)

# Fix a typo: atom_style is the correct keyword, not atom_sytle
write_lammps_data('structure.inp', atom, units='metal', atom_style='full')
