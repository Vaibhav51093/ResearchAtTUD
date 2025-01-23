from ovito import *
from ovito.io import import_file
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
import numpy as np
import sys
import matplotlib.pyplot as plt
import glob
import timeit

#path = "/nfshome/schulze/projects/20210215_NASICON/calculation/ID18_manual_twist/00-minimize/01-heat-523/02-diffusion"

#for submit script: for i in range(every_nth_frame, half_frame, every_nth_frame):

#Specifiy file path:
path = sys.argv[1]
#Reference frame
i = int(sys.argv[2])
every_animation_frame = int(sys.argv[3])

print(i)
# Chosing half-width of of grain boundary w in Angstrom
# Check with Ovito!
w = 10
#start = timeit.timeit()
###Import file


tot_list = glob.glob(f"{path}/structure.dump.*")
list = [file for file in tot_list if int(file.split(".dump.")[-1])%every_animation_frame == 0]  #Only every 10th dump-file loaded --> 500 of 5000 dump files 	
pipeline = import_file(list)
n_frames = pipeline.source.num_frames
Na_type_id = 5	

#timestep between frames in ps
#ÜBERPRÜFEN OB ZEITABSTÄNDE STIMMEN
dt = 0.002*100*(len(tot_list)/len(list))    # 2 ps pro Frame --> 500 Frames eingeladen
t=(i*dt)
#print(t, dt)

#Definitions for moving average
half_frame = pipeline.source.num_frames//2
#every_nth_frame = 1

# Set element

def setup_particle_types(frame, data):
    types = data.particles_.particle_types_
    types.type_by_id_(1).name = "Zr"
    types.type_by_id_(2).name = "Si"
    types.type_by_id_(3).name = "P"
    types.type_by_id_(4).name = "Na"
    types.type_by_id_(5).name = "O"
    types.type_by_id_(5).radius = 0.3
        
pipeline.modifiers.append(setup_particle_types)


# Set charges and mass

def set_charge_mass(frame: int, data: DataCollection):
    charge_dict = { "Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6, "Na":0.6, "Na":0.6}
    mass_dict = {"Zr": 91.224,
                        "Si": 28.086,
                        "P": 30.974,
                        "O": 15.999,
                        "Na": 22.99,
                        "Na": 22.99,
                        "Na": 22.99}

    charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
    mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
    ptypes = data.particles_.particle_types_
    for i in range(data.particles.count):
        charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
        mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]

pipeline.modifiers.append(set_charge_mass)


# Displacement vectors
displ = CalculateDisplacementsModifier(minimum_image_convention = False, use_frame_offset = True )
pipeline.modifiers.append(displ)


# Calculate 'Relative Displacement' of host matrix

def rel_displacement(frame, data, type_id_Na = 5):
        #with respect to center of mass of host matrix                                                                                                                                                                                        
        ptypes = data.particles['Particle Type']
        displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
        displ = data.particles['Displacement']
        masses = data.particles['Mass']
        #calculate center of mass displacement of host matrix                                                                                                                                                                                 
        total_mass_matrix = np.sum( masses[ (ptypes != type_id_Na) ] )
        sum = 0.0
        for i in range(len(displ_matrix)):
                sum +=displ_matrix[i] * masses[i]
        com_displ_host_matrix = sum/total_mass_matrix
        data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
        data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
import functools
pipeline.modifiers.append(functools.partial(rel_displacement, type_id_Na = 5))

###Custom python function that calculates and returns the msd_raw, netmsd_raw, msd_com, netmsd_com of atoms with a specific particle type id


########## neuMS
def define_gb_grain(frame: int, data: DataCollection, w = w, type_id_Na = 5):  ### "frame: int," benötigt?
        #Select GB: Choice has to be double-checked in Ovito!
        pos = data.particles.positions
        Selection = data.particles_.create_property("Selection")
        Selection[(pos[:,2] < data.cell[2][3] + w) ] = 1
        Selection[(pos[:,2] > data.cell[2][3] + data.cell[2][2]/2. - w) & (pos[:,2] < data.cell[2][3] + data.cell[2][2]/2. + w) ] = 1
        Selection[(pos[:,2] > data.cell[2][3] + data.cell[2][2] - w) ] = 1

        data.apply(FreezePropertyModifier(source_property = 'Selection', 
                                destination_property = 'Select0',
                                freeze_at = 0))

        ptype = data.particles.particle_type
        #only Na atoms
        Selection = data.particles["Select0"]
        #print(type(Selection))

pipeline.modifiers.append(define_gb_grain)


########## neuMS
def Calculate_MSD_gb_raw(data, type_id_Na):
        ptypes = data.particles["Particle Type"]
        Selection = data.particles["Select0"] #####neuMS
        gb_displ = data.particles['Displacement'][(ptypes == type_id_Na) & (Selection)]
        return np.mean(gb_displ**2, axis=0)

def Calculate_netMSD_gb_raw(data, type_id_Na):
        ptypes = data.particles["Particle Type"]
        Selection = data.particles["Select0"] #####neuMS
        gb_displ = data.particles['Displacement'][(ptypes == type_id_Na) & (Selection)]
        return np.mean(gb_displ, axis = 0)**2

def Calculate_MSD_gb_com(data, type_id_Na):
        ptypes = data.particles["Particle Type"]
        Selection = data.particles["Select0"] #####neuMS
        gb_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na) & (Selection)]
        return np.mean(gb_displ_rel**2, axis = 0)

def Calculate_netMSD_gb_com(data, type_id_Na):
        ptypes = data.particles["Particle Type"]
        Selection = data.particles["Select0"] #####neuMS
        gb_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na) & (Selection)]
        return np.mean(gb_displ_rel, axis = 0)**2


def Calculate_MSD_grain_raw(data, type_id_Na):
        ptypes = data.particles["Particle Type"]
        Selection = data.particles["Select0"] #####neuMS
        grain_displ = data.particles['Displacement'][(ptypes == type_id_Na) & (Selection != 1)]
        return np.mean(grain_displ**2, axis=0)

def Calculate_netMSD_grain_raw(data, type_id_Na):
        ptypes = data.particles["Particle Type"]
        Selection = data.particles["Select0"] #####neuMS
        grain_displ = data.particles['Displacement'][(ptypes == type_id_Na) & (Selection != 1)]
        return np.mean(grain_displ, axis = 0)**2

def Calculate_MSD_grain_com(data, type_id_Na):
        ptypes = data.particles["Particle Type"]
        Selection = data.particles["Select0"] #####neuMS
        grain_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na) & (Selection != 1)]
        return np.mean(grain_displ_rel**2, axis = 0)

def Calculate_netMSD_grain_com(data, type_id_Na):
        ptypes = data.particles["Particle Type"]
        Selection = data.particles["Select0"] #####neuMS
        grain_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na) & (Selection != 1)]
        return np.mean(grain_displ_rel, axis = 0)**2

# Smoothing
Na_MSD_gb = [np.array((0.0, 0.0, 0.0))]
net_Na_MSD_gb = [np.array((0.0, 0.0, 0.0))]
Na_MSD_com_gb = [np.array((0.0, 0.0, 0.0))]
net_Na_MSD_com_gb = [np.array((0.0, 0.0, 0.0))]
Na_MSD_grain = [np.array((0.0, 0.0, 0.0))]
net_Na_MSD_grain = [np.array((0.0, 0.0, 0.0))]
Na_MSD_com_grain = [np.array((0.0, 0.0, 0.0))]
net_Na_MSD_com_grain = [np.array((0.0, 0.0, 0.0))]



#print(f"0.000 \t {Na_MSD[-1][0]:.6e} \t {Na_MSD[-1][1]:.6e} \t {Na_MSD[-1][2]:.6e} \t {net_Na_MSD[-1][0]:.6e} \t {net_Na_MSD[-1][1]:.6e} \t {net_Na_MSD[-1][2]:.6e} \t {Na_MSD_com[-1][0]:.6e} \t {Na_MSD_com[-1][1]:.6e} \t {Na_MSD_com[-1][2]:.6e} \t {net_Na_MSD_com[-1][0]:.6e} \t {net_Na_MSD_com[-1][1]:.6e} \t {net_Na_MSD_com[-1][2]:.6e}")


displ.frame_offset = -i
Na_msd_gb = np.array((0.0, 0.0, 0.0))
net_Na_msd_gb = np.array((0.0, 0.0, 0.0))
Na_msd_com_gb = np.array((0.0, 0.0, 0.0))
net_Na_msd_com_gb = np.array((0.0, 0.0, 0.0))
Na_msd_grain = np.array((0.0, 0.0, 0.0))
net_Na_msd_grain = np.array((0.0, 0.0, 0.0))
Na_msd_com_grain = np.array((0.0, 0.0, 0.0))
net_Na_msd_com_grain = np.array((0.0, 0.0, 0.0))

for frame in range(i, n_frames): 
        data = pipeline.compute(frame)
        Na_msd_gb += Calculate_MSD_gb_raw(data, Na_type_id)
        net_Na_msd_gb += Calculate_netMSD_gb_raw(data, Na_type_id)
        Na_msd_com_gb += Calculate_MSD_gb_com(data, Na_type_id)
        net_Na_msd_com_gb += Calculate_netMSD_gb_com(data, Na_type_id)
        Na_msd_grain += Calculate_MSD_grain_raw(data, Na_type_id)
        net_Na_msd_grain += Calculate_netMSD_grain_raw(data, Na_type_id)
        Na_msd_com_grain += Calculate_MSD_grain_com(data, Na_type_id)
        net_Na_msd_com_grain += Calculate_netMSD_grain_com(data, Na_type_id)


Na_msd_gb/=(n_frames - i)
net_Na_msd_gb /=(n_frames -i )
Na_msd_com_gb/=(n_frames - i)
net_Na_msd_com_gb /=(n_frames -i )
Na_msd_grain/=(n_frames - i)
net_Na_msd_grain /=(n_frames -i )
Na_msd_com_grain/=(n_frames - i)
net_Na_msd_com_grain /=(n_frames -i )


#Export results to txt file

np.savetxt( f"{path}/03-smooth/diffusivity-attributes-rolling_mean_t{t:.0f}.txt", np.column_stack( (t, Na_msd_gb[0], Na_msd_gb[1], Na_msd_gb[2], net_Na_msd_gb[0], net_Na_msd_gb[1], net_Na_msd_gb[2], Na_msd_com_gb[0], Na_msd_com_gb[1], Na_msd_com_gb[2], net_Na_msd_com_gb[0], net_Na_msd_com_gb[1], net_Na_msd_com_gb[2], Na_msd_grain[0], Na_msd_grain[1], Na_msd_grain[2], net_Na_msd_grain[0], net_Na_msd_grain[1], net_Na_msd_grain[2], Na_msd_com_grain[0], Na_msd_com_grain[1], Na_msd_com_grain[2], net_Na_msd_com_grain[0], net_Na_msd_com_grain[1], net_Na_msd_com_grain[2] )), delimiter = " ")#, header = "Timestep(ps),MSDx_gb_raw(A^2),MSDy_gb_raw(A^2),MSDz_gb_raw(A^2),netMSDx_gb_raw(A^2),netMSDy_gb_raw(A^2),netMSDz_gb_raw(A^2),MSDx_gb_com(A^2),MSDy_gb_com(A^2),MSDz_gb_com(A^2),netMSDx_gb_com(A^2),netMSDy_gb_com(A^2),netMSDz_gb_com(A^2),MSDx_grain_raw(A^2),MSDy_grain_raw(A^2),MSDz_grain_raw(A^2),netMSDx_grain_raw(A^2),netMSDy_grain_raw(A^2),netMSDz_grain_raw(A^2),MSDx_grain_com(A^2),MSDy_grain_com(A^2),MSDz_grain_com(A^2),netMSDx_grain_com(A^2),netMSDy_grain_com(A^2),netMSDz_grain_com(A^2)")

