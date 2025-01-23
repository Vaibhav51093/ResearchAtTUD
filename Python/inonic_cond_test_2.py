import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import cumtrapz
import copy
import sys
import os
from ase.io import read
from fileinput import filename
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
import numpy as np 
from numba import jit, cuda

pipeline = import_file("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass/boilot_623K_2_05_cryt_523_k_ps_new_hdf5/boilot_623K_2_05_cryt_523_k_ps_new/dump_nvt_prod.out")

# Setting up the atom type 
pipeline.add_to_scene()
def setup_particle_types(frame, data):
    types = data.particles_.particle_types_
    types.type_by_id_(1).name = "O"
    types.type_by_id_(2).name = "Na"
    types.type_by_id_(3).name = "Zr"
    types.type_by_id_(4).name = "Si"
    types.type_by_id_(5).name = "P"
pipeline.modifiers.append(setup_particle_types)

def modify1(frame: int, data: DataCollection):

    charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
    mass_dict = {"Zr": 91.224,
                    "Si": 28.086,
                    "P": 30.974,
                    "O": 15.999,
                    "Na": 22.99}
    charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
    mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
    ptypes = data.particles_.particle_types_
    for i in range(data.particles.count):
        charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
        mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
pipeline.modifiers.append(modify1)

#def modify2(frame, data, type_id_Na = 2):
#    # with respect to center of mass of host matrix
#    ptypes = data.particles['Particle Type']
#    velocities_matrix = data.particles['Velocity'][ (ptypes != type_id_Na) ]
#    velocities = data.particles['Velocity']
#    masses = data.particles['Mass'][(ptypes != type_id_Na)]
#        
#    # calculate center of mass velocity of host matrix
#    total_mass_matrix = np.sum( masses)
#    sum = 0.0
#    for i in range(len(velocities_matrix)):
#        sum += velocities_matrix[i] * masses[i]
#    com_velocity_host_matrix = sum/total_mass_matrix
#    data.attributes['com_velocity_host_matrix'] = com_velocity_host_matrix
#    data.particles_.create_property('Relative Velocity', data = velocities - com_velocity_host_matrix)
#
#import functools
#pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))


#
def modify2(frame, data, type_id_Na = 2):
    # with respect to center of mass of host matrix
    ptypes = data.particles['Particle Type']
    velocities_matrix = data.particles['Velocity'][ (ptypes != type_id_Na) ]
    velocities = data.particles['Velocity']
    masses = data.particles['Mass'][(ptypes != type_id_Na)]
        
    # calculate center of mass velocity of host matrix
    total_mass_matrix = np.sum( masses)
    sum = 0.0
    for i in range(len(velocities_matrix)):
        sum += velocities_matrix[i] * masses[i]
    com_velocity_host_matrix = sum/total_mass_matrix
    data.attributes['com_velocity_host_matrix'] = com_velocity_host_matrix
    data.particles_.create_property('Relative Velocity', data = velocities - com_velocity_host_matrix)

import functools
pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

freeze_modifier = FreezePropertyModifier(source_property = 'Relative Velocity', destination_property = 'Velocity0', freeze_at = 0)
pipeline.modifiers.append(freeze_modifier)

# Define the function to calculate the ionic conductivity using green kubo method

def calculate_VAC(frame: int, data: DataCollection):
        ptypes = pipeline.particles.particle_type
        velocities_t0 = pipeline.particles["Velocity0"][(ptypes == 2)]
        velocities = pipeline.particles["Relative Velocity"][(ptypes == 2)]
        #sum_x = 0.0
        #sum_y = 0.0
        #sum_z = 0.0
        #for i in range(len(velocities)):
        #        sum_x += (velocities[i][0] * velocities_t0[i][0])
        #        sum_y += (velocities[i][1] * velocities_t0[i][1])
        #        sum_z += (velocities[i][2] * velocities_t0[i][2])
        sum_x/=len(velocities)
        sum_y/=len(velocities)
        sum_z/=len(velocities)
        #print(sum_x, sum_y, sum_z)                                                                                                                                                                                                    
        output.attributes["velocity-autocorr_x"] = sum_x
        output.attributes["velocity-autocorr_y"] = sum_y
        output.attributes["velocity-autocorr_z"] = sum_z
pipeline.modifiers.append(PythonScriptModifier(function = calculate_VAC))

mylist = []
for_export = []
interval_length = 500

# function optimized to run on gpu 
#@jit(target_backend='cuda')

for j in range(0, pipeline.source.num_frames-interval_length, int(interval_length/10)):
        print(j, j+interval_length, pipeline.source.num_frames)
        freeze_modifier.freeze_at = j
        for frame in range(j, j+interval_length):
                results = pipeline.compute(frame)
                for_export.append( [ frame, results.attributes["Timestep"], results.attributes["velocity-autocorr_x"],  results.attributes["velocity-autocorr_y"], results.attributes["velocity-autocorr_z"]]  )
        mylist.append(for_export)
        for_export = []

print(np.mean(mylist, axis = 0))
np.savetxt( "VAC_results_rolling_mean_2.txt", np.mean(mylist, axis = 0))

