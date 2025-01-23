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

# Run over several measurements

input_files = ["/nfshome/deshmukh/vaibhav/frank_jÜlich/lmp_6_1_3_2/simu_2/simu_2/523",
               "/nfshome/deshmukh/vaibhav/frank_jÜlich/lmp_6_1_3_2/simu_2/simu_2/573",
               "/nfshome/deshmukh/vaibhav/frank_jÜlich/lmp_6_1_3_2/simu_2/simu_2/623",
               "/nfshome/deshmukh/vaibhav/frank_jÜlich/lmp_6_1_3_2/simu_2/simu_2/673",
               "/nfshome/deshmukh/vaibhav/frank_jÜlich/lmp_6_1_3_2/simu_2/simu_2/723",
               "/nfshome/deshmukh/vaibhav/frank_jÜlich/lmp_6_1_3_2/simu_2/simu_2/773",]

# Chosing half-width of of grain boundary w in Angstrom
#w = 10 # Check with Ovito!

for element in input_files:
    # Data import
    infile = element
    
    os.chdir = infile
    pipeline = import_file(f"{infile}/dump_nvt_prod.out")
    print(infile)

    POTIM = 1                  # Time step in fs
    NBLOCK = 100            # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
    # MSD Evaluation
    Stepsize = 10
    # Set element

    def setup_particle_types(frame, data):
        types = data.particles_.particle_types_
        types.type_by_id_(1).name = "O"
        types.type_by_id_(2).name = "Na"
        types.type_by_id_(4).name = "Si"
        types.type_by_id_(5).name = "P"
        types.type_by_id_(5).radius = 0.3
        
    pipeline.modifiers.append(setup_particle_types)


    # Set charges and mass

    def modify1(frame: int, data: DataCollection):

        charge_dict = { "Si": 2.4, "P": 3.0, "O":-1.2, "Na":0.6, "Na":0.6, "Na":0.6}
        mass_dict = { "Si": 28.086,
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

    pipeline.modifiers.append(modify1)


    # Displacement vectors

    pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

    # Minimum image convention is not used, because the system is periodic in all directions and also unwrapped 

    # Calculate 'Relative Displacement' of host matrix

    def modify2(frame, data, type_id_Na = 2):
            
        #with respect to center of mass of host matrix                                                                                                                                                                                        
        ptypes = data.particles['Particle Type']
        displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
        displ = data.particles['Displacement']
        masses = data.particles['Mass'][(ptypes != type_id_Na)]
            
        #calculate center of mass displacement of host matrix                                                                                                                                                                                 
        total_mass_matrix = np.sum( masses)
        sum = 0.0
        for i in range(len(displ_matrix)):
                sum +=displ_matrix[i] * masses[i]
        com_displ_host_matrix = sum/total_mass_matrix
        data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
        data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
    
    import functools
    pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

    # Calculate MSD, netMSD for x-,y-,z-direction

    def modify3(frame: int, data: DataCollection, type_id_Na = 2):

        ptype = data.particles.particle_type

        time = frame*NBLOCK*POTIM/1000

        pos = data.particles['Position'][(ptype == type_id_Na)]
        #print(len(pos[:,0]))
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        bulk_displ = data.particles['Displacement'][(ptype == type_id_Na)]
        #print(bulk_displ_rel)
        #rel_x, rel_y, rel_z = np.mean(bulk_displ_rel, axis = 0)
        
        msd_x, msd_y, msd_z = np.mean(bulk_displ_rel**2, axis = 0)
        data.attributes["MSDx_bulk_com"] = np.mean(bulk_displ_rel[:, 0] * bulk_displ_rel[:, 0], axis=0)  # com = center of mass /displacements are corrected for com shift of host matrix
        data.attributes["MSDy_bulk_com"] = np.mean(bulk_displ_rel[:, 1] * bulk_displ_rel[:, 1], axis=0)  # MSD = tracer diffucivity
        data.attributes["MSDz_bulk_com"] = np.mean(bulk_displ_rel[:, 2] * bulk_displ_rel[:, 2], axis=0)
        data.attributes["MSDxy_bulk_com"] = np.mean(bulk_displ_rel[:, 0] * bulk_displ_rel[:, 1], axis=0) #np.mean(cd_xy**2) #np.mean(np.array(cd_xy)**2, axis=1) 
        data.attributes["MSDxz_bulk_com"] = np.mean(bulk_displ_rel[:, 0] * bulk_displ_rel[:, 2], axis=0) #np.mean(cd_xz**2) #np.mean(np.array(cd_xz)**2, axis=1) 
        data.attributes["MSDyz_bulk_com"] = np.mean(bulk_displ_rel[:, 1] * bulk_displ_rel[:, 2], axis=0)
        
        data.attributes["MSDav_bulk_com"] = msd_x+msd_y+msd_z
        
        # The follwoing code implement cros_displacement algo to get offdiagonal compoenets of MSD 
        #sum_of_pos_xy = np.mean((bulk_displ_rel[:,0] + bulk_displ_rel[:,1])**2)
        #sum_of_pos_xz = np.mean((bulk_displ_rel[:,0] + bulk_displ_rel[:,2])**2)
        #sum_of_pos_yz = np.mean((bulk_displ_rel[:,1] + bulk_displ_rel[:,2])**2)
        
        #cd_xy = cd[0][1]
        #cd_xz = cd[0][2]
        #cd_yz = cd[1][2]
        
        #print(cd_xy)
        
        #data.attributes["MSDxy_bulk_com"] = 0.5*(sum_of_pos_xy - msd_x - msd_y) #np.mean(cd_xy**2) #np.mean(np.array(cd_xy)**2, axis=1) 
        #data.attributes["MSDxz_bulk_com"] = 0.5*(sum_of_pos_xz - msd_x - msd_z) #np.mean(cd_xz**2) #np.mean(np.array(cd_xz)**2, axis=1) 
        #data.attributes["MSDyz_bulk_com"] = 0.5*(sum_of_pos_yz - msd_y - msd_z) #np.mean(cd_yz**2) #np.mean(np.array(cd_yz)**2, axis=1)

        net_x, net_y, net_z = np.sum(bulk_displ_rel, axis = 0)
        data.attributes["netMSDx_bulk_com"] = net_x**2 / len(bulk_displ_rel) #netMSD = ionic diffusivity
        data.attributes["netMSDy_bulk_com"] = net_y**2 / len(bulk_displ_rel)
        data.attributes["netMSDz_bulk_com"] = net_z**2 / len(bulk_displ_rel)
        data.attributes["netMSDav_bulk_com"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ_rel)
        
        # The follwoing code implement cros_displacement algo to get offdiagonal compoenets of MSD 
        sum_of_pos_xy_n = ((np.sum(bulk_displ_rel[:,0]) + np.sum(bulk_displ_rel[:,1]))**2) / len(bulk_displ_rel)
        sum_of_pos_xz_n = ((np.sum(bulk_displ_rel[:,0]) + np.sum(bulk_displ_rel[:,2]))**2) / len(bulk_displ_rel)
        sum_of_pos_yz_n = ((np.sum(bulk_displ_rel[:,1]) + np.sum(bulk_displ_rel[:,2]))**2) / len(bulk_displ_rel)
        
        data.attributes["netMSDxy_bulk_com"] = 0.5*(sum_of_pos_xy_n - net_x**2/len(bulk_displ_rel) - net_y**2/len(bulk_displ_rel)) #np.mean(cd_xy**2) #np.mean(np.array(cd_xy)**2, axis=1) 
        data.attributes["netMSDxz_bulk_com"] = 0.5*(sum_of_pos_xz_n - net_x**2/len(bulk_displ_rel) - net_z**2/len(bulk_displ_rel)) #np.mean(cd_xz**2) #np.mean(np.array(cd_xz)**2, axis=1) 
        data.attributes["netMSDyz_bulk_com"] = 0.5*(sum_of_pos_yz_n - net_y**2/len(bulk_displ_rel) - net_z**2/len(bulk_displ_rel))
        
        data.attributes["Timestep"] = time
     
    pipeline.modifiers.append(modify3)
    
    export_file(pipeline, f"{infile}/msd_data_3_1_3.txt", "txt/attr", columns=["Timestep", "MSDx_bulk_com", "MSDy_bulk_com", "MSDz_bulk_com", "MSDav_bulk_com","MSDxy_bulk_com","MSDxz_bulk_com","MSDyz_bulk_com","netMSDx_bulk_com", "netMSDy_bulk_com", "netMSDz_bulk_com","netMSDav_bulk_com", "netMSDxy_bulk_com", "netMSDxz_bulk_com", "netMSDyz_bulk_com"],multiple_frames=True)

    
