import os
import numpy as np
from ovito.io import import_file, export_file
from ovito.pipeline import Pipeline
from ovito.data import DataCollection

# Data import
infile = 'dump_nvt_523_Na6Si3P4O19.out'  # Replace with your directory name
os.chdir(infile)
pipeline = import_file(f"dump_nvt_523_Na6Si3P4O19_2.out")
print(infile)

# Simulation settings
POTIM = 1                  # Time step in femtoseconds
NBLOCK = 100               # Steps per output frame in the simulation
Stepsize = 10              # Export interval in picoseconds

# Setup particle types and properties
def setup_particle_types(frame, data):
    types = data.particles.particle_types
    types.type_by_id(1).name = "O"
    types.type_by_id(2).name = "Na"
    types.type_by_id(4).name = "Si"
    types.type_by_id(5).name = "P"
    types.type_by_id(5).radius = 0.3

pipeline.modifiers.append(setup_particle_types)

def modify1(frame, data):
    charge_dict = {"Si": 2.4, "P": 3.0, "O": -1.2, "Na": 0.6}
    mass_dict = {"Si": 28.086, "P": 30.974, "O": 15.999, "Na": 22.99}
    charge = data.particles_.create_property("Charge", data=np.ones(data.particles.count))
    mass = data.particles_.create_property("Mass", data=np.ones(data.particles.count))
    ptypes = data.particles.particle_types
    for i in range(data.particles.count):
        ptype_name = ptypes.type_by_id(ptypes[i]).name
        charge[i] = charge_dict[ptype_name]
        mass[i] = mass_dict[ptype_name]

pipeline.modifiers.append(modify1)

# Determine the output frequency in timesteps
output_frequency_steps = Stepsize * (10000 / POTIM)  # Convert ps to steps

# Export new file every specified interval
def export_conditionally(frame, data):
    # Check if the current frame should trigger an export
    if frame % output_frequency_steps == 0:
        timestep = data.attributes['Timestep']
        export_file(data, f'output_frame_{timestep}.dump', 'lammps/dump',
                    columns=["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])

pipeline.modifiers.append(export_conditionally)

# Process all frames
pipeline.compute()

