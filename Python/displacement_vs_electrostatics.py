# ----------------------------------------------------------------------
# This is the official script of Vaibhav Arun Deshmukh
# Script is written in Python 3.7.4
# It requires Ovito 3.0.0-dev454 or higher
# ----------------------------------------------------------------------
# Purpose: This script is used to calculate the displacement of ions, 
#          and the electrostatic potential at the same time, plus the
#          local environments of the ions (10 Angstroms)
# I will try to save data in csv format, later on I will try to plot 
# -----------------------------------------------------------------------
from fileinput import filename
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
import numpy as np 
import numpy as np
from scipy import integrate
from scipy.constants import codata
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
plt.style.use(['vaibhz-sci','no-latex','high-vis','vibrant'])
# -------------------------------------------------------------------------
# Local analysis and electrostatics loccal analysis, Enjoy!
# -------------------------------------------------------------------------

print('Starting the script')

pipeline_1 = import_file("/nfshome/deshmukh/vaibhav/shanka_cluster/deshmukh2/nasi_2_random/623/dump_nvt_prod_3.out", multiple_frames = True)
pipeline_2 = import_file("/nfshome/deshmukh/vaibhav/shanka_cluster/deshmukh2/nasi_2_random/623/dump_nvt_prod_3.out", multiple_frames = True)

# Basic info
n_frames = pipeline_1.source.num_frames    # Number of frames in the trajectory
Na_type_id = 2                           # Type id of Na
POTIM = 2.0                              # Time step in fs
NBLOCK = 200                             # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
Stepsize = 10                            # How many steps to skip in the trajectory

# Add particle types to the scene
pipeline_1.add_to_scene()
def setup_particle_types(frame, data):
    types = data.particles_.particle_types_
    types.type_by_id_(1).name = "O"
    types.type_by_id_(2).name = "Na"
    types.type_by_id_(3).name = "Zr"
    types.type_by_id_(4).name = "Si"
    types.type_by_id_(5).name = "P"
pipeline_1.modifiers.append(setup_particle_types)

# Add mass and charge to the particles
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
        
pipeline_1.modifiers.append(modify1)

data = pipeline_1.compute()                                                               # Compute the pipeline
ptypes = data.particles.particle_types  

# Move the modifier outside the loop if possible
#pipeline.modifiers.append(WrapPeriodicImagesModifier())
#mod = CalculateDisplacementsModifier(use_frame_offset = True, frame_offset = -25)
mod = CalculateDisplacementsModifier(use_frame_offset = False) #, frame_offset = -2)
pipeline_1.modifiers.append(mod)
# Displacement vectors:
# Initialize lists to store displacements

avg_disp_x_0 = []
avg_disp_y_0 = []
avg_disp_z_0 = []
avg_total_displacement_0 = []

frame_queue_0 = []

for j in range(0, 25000, 1):
    data = pipeline_1.compute(j)
    na_index = data.particles["Particle Identifier"][(data.particles.particle_type == 2)]
    
    disp_x = []
    disp_y = []
    disp_z = []
    total_displacement = []
    
    #frame_queue_0 = []
    
    for k in range(len(na_index)):
        d = data.particles['Displacement'][na_index[k]]
        
        disp_x.append(d[0])
        disp_y.append(d[1])
        disp_z.append(d[2])
        total_displacement.append(np.sqrt(d[0]**2 + d[1]**2 + d[2]**2))#/3.4778) # We are calculating hoping factor 
    
    # Add the current frame's data to the frame queue
    frame_queue_0.append({
        'disp_x': disp_x,
        'disp_y': disp_y,
        'disp_z': disp_z,
        'total_displacement': total_displacement
    })
    
    # If the frame queue has more than 25 frames, remove the oldest frame
    if len(frame_queue_0) > 25:
        removed_frame = frame_queue_0.pop(0)
        
    # Calculate the average displacement over the last 25 frames
    if len(frame_queue_0)==25:
        #avg_si_count_max = [sum(x) / 25 for x in zip(*[frame['si_count_max_1'] for frame in frame_queue])]
        avg_disp_x = [sum(x) / 25 for x in zip(*[frame['disp_x'] for frame in frame_queue_0])]
        avg_disp_y = [sum(x) / 25 for x in zip(*[frame['disp_y'] for frame in frame_queue_0])]
        avg_disp_z = [sum(x) / 25 for x in zip(*[frame['disp_z'] for frame in frame_queue_0])]
        avg_total_displacement = [sum(x) / 25 for x in zip(*[frame['total_displacement'] for frame in frame_queue_0])]

        # append the average to the list
        avg_disp_x_0.append(avg_disp_x)
        avg_disp_y_0.append(avg_disp_y)
        avg_disp_z_0.append(avg_disp_z)
        avg_total_displacement_0.append(avg_total_displacement)

d_x = np.array(avg_disp_x_0).T
d_y = np.array(avg_disp_y_0).T
d_z = np.array(avg_disp_z_0).T
total_d = np.array(avg_total_displacement_0).T

cutoff = 5

# Initialize lists to store data for the last 25 frames
index_si_max_2 = []
index_p_max_2 = []
index_na_max_2 = []
index_zr_max_2 = []

si_count_max_2 = []
p_count_max_2 = []
na_count_max_2 = []
zr_count_max_2 = []

# Initialize a queue to keep track of the last 25 frames
frame_queue = []

for i in range(0, 25000, 1):
    data = pipeline_2.compute(i)
    na_index = data.particles["Particle Identifier"][(data.particles.particle_type == 2)]
    partcle_type_max = data.particles.particle_types[na_index]
    finder = CutoffNeighborFinder(cutoff, data)

    positions = data.particles.positions
    neighbors = finder.find(na_index)
    
    si_count_max_1 = []
    p_count_max_1 = []
    na_count_max_1 = []
    zr_count_max_1 = []

    index_si_max_1 = []
    index_p_max_1 = []
    index_na_max_1 = []
    index_zr_max_1 = []
    
    for k in range(len(na_index)):
        si_count_max = 0
        p_count_max = 0
        na_count_max = 0
        zr_count_max = 0

        index_si_max = []
        index_p_max = []
        index_na_max = []
        index_zr_max = []
              
        for neigh in finder.find(na_index[k]):
            if data.particles.particle_types[neigh.index] == 4:              
                index_si_max.append(neigh.index)
                si_count_max += 1
            elif data.particles.particle_types[neigh.index] == 5:
                index_p_max.append(neigh.index)
                p_count_max += 1
            elif data.particles.particle_types[neigh.index] == 2:
                index_na_max.append(neigh.index)
                na_count_max += 1
            elif data.particles.particle_types[neigh.index] == 3:
                index_zr_max.append(neigh.index)
                zr_count_max += 1

        si_count_max_1.append(si_count_max)
        p_count_max_1.append(p_count_max)
        na_count_max_1.append(na_count_max)
        zr_count_max_1.append(zr_count_max)

        index_si_max_1.append(index_si_max)
        index_p_max_1.append(index_p_max)
        index_na_max_1.append(index_na_max)
        index_zr_max_1.append(index_zr_max)

    # Add the current frame's data to the frame queue
    frame_queue.append({
        'si_count_max_1': si_count_max_1,
        'p_count_max_1': p_count_max_1,
        'na_count_max_1': na_count_max_1,
        'zr_count_max_1': zr_count_max_1,
        'index_si_max_1': index_si_max_1,
        'index_p_max_1': index_p_max_1,
        'index_na_max_1': index_na_max_1,
        'index_zr_max_1': index_zr_max_1
    })

    # If the frame queue has more than 25 frames, remove the oldest frame
    if len(frame_queue) > 25:
        removed_frame = frame_queue.pop(0)

    # Calculate the average over the last 25 frames
    if len(frame_queue) == 25:
        avg_si_count_max = [sum(x) / 25 for x in zip(*[frame['si_count_max_1'] for frame in frame_queue])]
        avg_p_count_max = [sum(x) / 25 for x in zip(*[frame['p_count_max_1'] for frame in frame_queue])]
        avg_na_count_max = [sum(x) / 25 for x in zip(*[frame['na_count_max_1'] for frame in frame_queue])]
        avg_zr_count_max = [sum(x) / 25 for x in zip(*[frame['zr_count_max_1'] for frame in frame_queue])]

        si_count_max_2.append(avg_si_count_max)
        p_count_max_2.append(avg_p_count_max)
        na_count_max_2.append(avg_na_count_max)
        zr_count_max_2.append(avg_zr_count_max)
        
# x, potential

charge_density_x = []
potential_x_list = []
data_x = []

frame_queue_1 = []

for j in range(0, 25000, 1):
    charge_density_x_1 = []
    potential_x_list_1 = []
    data_x_1 = []

    for k in range(len(na_index)):
        # Clear modifiers specific to this particle
        pipeline_1.modifiers.clear()
        
        # Add particle-specific modifier
        pipeline_1.modifiers.append(ExpressionSelectionModifier(expression='ParticleIdentifier == {}'.format(na_index[k])))
        pipeline_1.modifiers.append(ExpandSelectionModifier(cutoff=5.0))
        pipeline_1.modifiers.append(InvertSelectionModifier())
        pipeline_1.modifiers.append(DeleteSelectedModifier())
        pipeline_1.modifiers.append(ComputePropertyModifier(expressions = ('1',), output_property = 'Charge'))
        pipeline_1.modifiers.append(SpatialBinningModifier(property='Charge', reduction_operation=SpatialBinningModifier.Operation.SumVol, direction=SpatialBinningModifier.Direction.X, bin_count=(200, 80, 46)))

        # Compute data using the modified pipeline_1
        data = pipeline_1.compute(j)
        cell_x = data.cell[:,0][0]
        ch_density_x = data.tables['binning']['Charge']

        rho_x = ch_density_x
        e_field_x = 14.3997584 * integrate.cumtrapz(rho_x, np.linspace(0, cell_x, 200), initial=0)
        e_field_x = e_field_x - np.mean(e_field_x)
        potential_x = -integrate.cumtrapz(e_field_x, np.linspace(0, cell_x, 200))

        potential_x_scalar = np.array([potential_x[0]])
        #potential_x_combined = np.concatenate(potential_x_scalar)
        potential_x_list_1.append(potential_x_scalar)
        data_x_1.append(np.linspace(0, cell_x, 200))
        
    frame_queue_1.append({
        'potential_x_list_1': potential_x_list_1,
        'charge_density_x_1': charge_density_x_1,
        'data_x_1': data_x_1
    })
    
    # If the frame queue has more than 25 frames, remove the oldest frame
    if len(frame_queue_1) > 25:
        removed_frame = frame_queue_1.pop(0)
        
    # Calculate the average displacement over the last 25 frames
    if len(frame_queue_1)==25:
        
        avg_pot_x = [sum(x) / 25 for x in zip(*[frame['potential_x_list_1'] for frame in frame_queue_1])]
        avg_charge_density_x = [sum(x) / 25 for x in zip(*[frame['charge_density_x_1'] for frame in frame_queue_1])]
        avg_data_x = [sum(x) / 25 for x in zip(*[frame['data_x_1'] for frame in frame_queue_1])]

        charge_density_x.append(avg_charge_density_x)
        potential_x_list.append(avg_pot_x)
        data_x.append(avg_data_x)
        
# y, potential

charge_density_y = []
potential_y_list = []
data_y = []

frame_queue_2 = []

for j in range(0, 25000, 1):
    
    charge_density_y_1 = []
    potential_y_list_1 = []
    data_y_1 = []

    for k in range(len(na_index)):
        # Clear modifiers specific to this particle
        pipeline_1.modifiers.clear()
        
        # Add particle-specific modifier
        pipeline_1.modifiers.append(ExpressionSelectionModifier(expression='ParticleIdentifier == {}'.format(na_index[k])))
        pipeline_1.modifiers.append(ExpandSelectionModifier(cutoff=5.0))
        pipeline_1.modifiers.append(InvertSelectionModifier())
        pipeline_1.modifiers.append(DeleteSelectedModifier())
        pipeline_1.modifiers.append(ComputePropertyModifier(expressions = ('1',), output_property = 'Charge'))
        pipeline_1.modifiers.append(SpatialBinningModifier(property='Charge', reduction_operation=SpatialBinningModifier.Operation.SumVol, direction=SpatialBinningModifier.Direction.Y, bin_count=(200, 80, 46)))

        # Compute data using the modified pipeline_1
        data = pipeline_1.compute(j)
        cell_y = data.cell[:,1][1]
        ch_density_y = data.tables['binning']['Charge']

        rho_y = ch_density_y
        e_field_y = 14.3997584 * integrate.cumtrapz(rho_y, np.linspace(0, cell_y, 200), initial=0)
        e_field_y = e_field_y - np.mean(e_field_y)
        potential_y = -integrate.cumtrapz(e_field_y, np.linspace(0, cell_y, 200))

        potential_y_scalar = np.array([potential_y[0]])
        #potential_y_combined = np.concatenate(potential_y_scalar)
        potential_y_list_1.append(potential_y_scalar)
        data_y_1.append(np.linspace(0, cell_y, 200))
        
    frame_queue_2.append({
        'potential_y_list_1': potential_y_list_1,
        'charge_density_y_1': charge_density_y_1,
        'data_y_1': data_y_1
    })
    
    # If the frame queue has more than 25 frames, remove the oldest frame
    if len(frame_queue_2) > 25:
        removed_frame = frame_queue_2.pop(0)
        
    # Calculate the average displacement over the last 25 frames
    if len(frame_queue_2)==25:
        
        avg_pot_y = [sum(x) / 25 for x in zip(*[frame['potential_y_list_1'] for frame in frame_queue_2])]
        avg_charge_density_y = [sum(x) / 25 for x in zip(*[frame['charge_density_y_1'] for frame in frame_queue_2])]
        avg_data_y = [sum(x) / 25 for x in zip(*[frame['data_y_1'] for frame in frame_queue_2])]

        charge_density_y.append(avg_charge_density_y)
        potential_y_list.append(avg_pot_y)
        data_y.append(avg_data_y)
        

    charge_density_y.append(charge_density_y_1)
    potential_y_list.append(potential_y_list_1)
    data_y.append(data_y_1)
    
# z, potential

charge_density_z = []
potential_z_list = []
data_z = []

frame_queue_3 = []

for j in range(0, 25000, 1):
    
    charge_density_z_1 = []
    potential_z_list_1 = []
    data_z_1 = []

    for k in range(len(na_index)):
        # Clear modifiers specific to this particle
        pipeline_1.modifiers.clear()
        
        # Add particle-specific modifier
        pipeline_1.modifiers.append(ExpressionSelectionModifier(expression='ParticleIdentifier == {}'.format(na_index[k])))
        pipeline_1.modifiers.append(ExpandSelectionModifier(cutoff=5.0))
        pipeline_1.modifiers.append(InvertSelectionModifier())
        pipeline_1.modifiers.append(DeleteSelectedModifier())
        pipeline_1.modifiers.append(ComputePropertyModifier(expressions = ('1',), output_property = 'Charge'))
        pipeline_1.modifiers.append(SpatialBinningModifier(property='Charge', reduction_operation=SpatialBinningModifier.Operation.SumVol, direction=SpatialBinningModifier.Direction.Z, bin_count=(200, 80, 46)))

        # Compute data using the modified pipeline_1
        data = pipeline_1.compute(j)
        cell_z = data.cell[:,2][2]
        ch_density_z = data.tables['binning']['Charge']

        rho_z = ch_density_z
        e_field_z = 14.3997584 * integrate.cumtrapz(rho_z, np.linspace(0, cell_z, 200), initial=0)
        e_field_z = e_field_z - np.mean(e_field_z)
        potential_z = -integrate.cumtrapz(e_field_z, np.linspace(0, cell_z, 200))

        potential_z_scalar = np.array([potential_z[0]])
        #potential_z_combined = np.concatenate(potential_z_scalar)
        potential_z_list_1.append(potential_z_scalar)
        data_z_1.append(np.linspace(0, cell_z, 200))
        
    frame_queue_3.append({
        'potential_z_list_1': potential_z_list_1,
        'charge_density_z_1': charge_density_z_1,
        'data_z_1': data_z_1
    })
    
    # If the frame queue has more than 25 frames, remove the oldest frame
    if len(frame_queue_3) > 25:
        removed_frame = frame_queue_3.pop(0)
        
    # Calculate the average displacement over the last 25 frames
    if len(frame_queue_3)==25:
            
        avg_pot_z = [sum(x) / 25 for x in zip(*[frame['potential_z_list_1'] for frame in frame_queue_3])]
        avg_charge_density_z = [sum(x) / 25 for x in zip(*[frame['charge_density_z_1'] for frame in frame_queue_3])]
        avg_data_z = [sum(x) / 25 for x in zip(*[frame['data_z_1'] for frame in frame_queue_3])]
    
        charge_density_z.append(avg_charge_density_z)
        potential_z_list.append(avg_pot_z)
        data_z.append(avg_data_z)    
    
def moving_average_numpy(data, window):
    weights = np.repeat(1.0, window) / window
    convolved_data = np.convolve(data, weights, mode='valid')
    return convolved_data

# Example data (replace this with your actual data)

window_size = 5

x = d_x.T
y = d_y.T
z = d_z.T

na_index = data.particles["Particle Identifier"][(data.particles.particle_type == 2)]

# Combine indices and displacement magnitudes
mx = []
#mn = []

for i, j in zip(na_index, na_index):
    mx.append(np.where(na_index == i)[0][0])
    #mn.append(np.where(na_index == j)[0][0])

mx = np.array(mx)
#mn = np.array(mn)

dis = total_d.T

def running_average(si_min):
  """Calculates the running average of a list."""
  running_sum = 0
  count = 0
  running_averages = []
  for si_min_value in si_min:
    running_sum += si_min_value
    count += 1
    running_average = running_sum / count
    running_averages.append(running_average)
  return running_averages

import math

def running_average(si_min):
    """Calculates the running average of a list."""
    running_sum = 0
    count = 0
    running_averages = []
    for si_min_value in si_min:
        running_sum = math.fsum([running_sum, si_min_value])
        count += 1
        running_average = running_sum / count
        running_averages.append(running_average)
    return running_averages

combined_data_p = list(zip(p_count_max_2, list(total_displacement)))
combined_data_si = list(zip(si_count_max_2, list(total_displacement)))
combined_data_na = list(zip(na_count_max_2, list(total_displacement)))
combined_data_zr = list(zip(zr_count_max_2, list(total_displacement)))

potential_x_list_2 = np.array(potential_x_list).T
potential_x_list_3 = potential_x_list_2[0].T

potential_y_list_2 = np.array(potential_y_list).T
potential_y_list_3 = potential_y_list_2[0].T

potential_z_list_2 = np.array(potential_z_list).T
potential_z_list_3 = potential_z_list_2[0].T

all_combined_adata = list(zip(p_count_max_2, si_count_max_2, na_count_max_2, zr_count_max_2, x, y, z, dis , potential_x_list_3, potential_y_list_3, potential_z_list_3))

sorted_data_p = sorted(combined_data_p, key=lambda x: x[0])
sorted_data_si = sorted(combined_data_si, key=lambda x: x[0])
sorted_data_na = sorted(combined_data_na, key=lambda x: x[0])
sorted_data_zr = sorted(combined_data_zr, key=lambda x: x[0])

p_max = np.array(p_count_max_2)   #.T
si_max = np.array(si_count_max_2) #.T
na_max = np.array(na_count_max_2) #.T
zr_max = np.array(zr_count_max_2) #.T

import numpy as np

# Assuming all_combined_adata is a list of tuples (p_count_max_2, si_count_max_2, na_count_max_2, zr_count_max_2, x, y, z, total_displacement)

# Separate the data into individual lists
p_count_max_2, si_count_max_2, na_count_max_2, zr_count_max_2, x, y, z, total_displacement, pot_x, pot_y, pot_z = zip(*all_combined_adata)

# Convert p_count_max_2 to a numpy array for sorting
p_count_max_2 = np.array(p_count_max_2)

# Find unique values of p_count_max_2
unique_p_count_max_2 = np.unique(p_count_max_2)

# Initialize lists to store the averaged values
averaged_si_count_max_2 = []
averaged_na_count_max_2 = []
averaged_zr_count_max_2 = []
averaged_x = []
averaged_y = []
averaged_z = []
averaged_total_displacement = []
averaged_total_displacement_2 = []
pot_x_t = []
pot_x_t_2 = []
pot_y_t = []
pot_y_t_2 = []
pot_z_t = []
pot_z_t_2 = []

# Iterate through unique values of p_count_max_2 and average the corresponding data
for unique_p in unique_p_count_max_2:
    mask = (p_count_max_2 == unique_p)
    averaged_si_count_max_2.append(np.mean(np.array(si_count_max_2)[mask]))
    averaged_na_count_max_2.append(np.mean(np.array(na_count_max_2)[mask]))
    averaged_zr_count_max_2.append(np.mean(np.array(zr_count_max_2)[mask]))
    averaged_x.append(np.mean(np.array(x)[mask]))
    averaged_y.append(np.mean(np.array(y)[mask]))
    averaged_z.append(np.mean(np.array(z)[mask]))
    averaged_total_displacement.append(np.array(total_displacement)[mask])
    averaged_total_displacement_2.append(np.mean(np.array(total_displacement)[mask]))
    pot_x_t.append(np.array(pot_x)[mask])
    pot_x_t_2.append(np.mean(np.array(pot_x)[mask]))
    pot_y_t.append(np.array(pot_y)[mask])
    pot_y_t_2.append(np.mean(np.array(pot_y)[mask]))
    pot_z_t.append(np.array(pot_z)[mask])
    pot_z_t_2.append(np.mean(np.array(pot_z)[mask]))

# Combine the averaged values back into a list of tuples
averaged_all_combined_adata = list(zip(
    unique_p_count_max_2,
    averaged_si_count_max_2,
    averaged_na_count_max_2,
    averaged_zr_count_max_2,
    averaged_x,
    averaged_y,
    averaged_z,
    averaged_total_displacement,
    averaged_total_displacement_2,
    pot_x_t,
    pot_x_t_2,
    pot_y_t,
    pot_y_t_2,
    pot_z_t,
    pot_z_t_2
))

# Sort the averaged_all_combined_adata based on unique_p_count_max_2 (if needed)
sorted_averaged_all_combined_adata = sorted(averaged_all_combined_adata, key=lambda x: x[0])

# Print or use sorted_averaged_all_combined_adata as needed

unique_p_count_max_2, averaged_si_count_max_2, averaged_na_count_max_2, averaged_zr_count_max_2, x_, y_, z_, averaged_total_displacement,averaged_total_displacement_2, p_x, p_x_2, p_y, p_y_2, p_z, p_z_2= zip(*averaged_all_combined_adata)

# Create a separate plot for each column
fig, axs = plt.subplots(5, 1, figsize=(10,9), sharex=True)

# Plot each column against unique_p_count_max_2
axs[0].plot(unique_p_count_max_2, averaged_si_count_max_2, marker='o')
axs[0].set_ylabel('Si ')
axs[0].legend()
axs[0].grid()

axs[1].plot(unique_p_count_max_2, averaged_na_count_max_2, marker='o')
axs[1].set_ylabel('Na')
axs[1].legend()
axs[1].grid()

axs[2].plot(unique_p_count_max_2, averaged_zr_count_max_2, marker='o')
axs[2].set_ylabel('Zr')

axs[2].legend()

axs[2].grid()

axs[3].plot(unique_p_count_max_2, averaged_total_displacement_2, marker='o')
axs[3].set_ylabel('Displacement')
axs[3].legend()
axs[3].grid()

axs[4].plot(unique_p_count_max_2, p_x_2, marker='o')
axs[4].set_ylabel('Potential')
axs[4].legend()
axs[4].grid()

# Adjust subplot layout
plt.tight_layout()
plt.xlabel('P')

# Show the plots
plt.savefig('si_na_vs_p_count_max_2_random.png')
plt.close()
#plt.show()

p_tot = np.array(p_x) + np.array(p_y) + np.array(p_z)
p_tot_2 = np.array(p_x_2) + np.array(p_y_2) + np.array(p_z_2)

# Create a DataFrame
data = pd.DataFrame({
    'P': unique_p_count_max_2,
    'Si': averaged_si_count_max_2,
    'Na': averaged_na_count_max_2,
    'Zr': averaged_zr_count_max_2,
    'Displacement': averaged_total_displacement_2,
    'Potential': p_tot_2
})

# Create a scatter plot matrix with marginal histograms and correlations
sns.set(style="ticks")
g = sns.pairplot(data, diag_kind="kde", markers="o", height=2)

# Add correlation coefficients to the upper triangle of the matrix
corr_matrix = data.corr()
for i, j in zip(*np.triu_indices_from(g.axes, 1)):
    g.axes[i, j].annotate(
        f"corr = {corr_matrix.iloc[i, j]:.2f}",
        (0.5, 0.95),
        xycoords="axes fraction",
        ha="center",
        va="center",
        fontsize=10,
    )

plt.suptitle("Scatter Plot Matrix with Marginal Histograms and Correlations", y=1.02)
plt.savefig('si_na_vs_p_count_max_2_random_2.png')
plt.close()

print('Done')