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

input_files = ["/nfshome/deshmukh/vaibhav/schulze_project/gb/IDp03_aimsgbCK/00-minimize/01-heat-523/02-diffusion",
               "/nfshome/deshmukh/vaibhav/schulze_project/gb/IDp03_aimsgbCK/00-minimize/01-heat-573/02-diffusion",
               "/nfshome/deshmukh/vaibhav/schulze_project/gb/IDp03_aimsgbCK/00-minimize/01-heat-623/02-diffusion",
               "/nfshome/deshmukh/vaibhav/schulze_project/gb/IDp03_aimsgbCK/00-minimize/01-heat-673/02-diffusion",
               "/nfshome/deshmukh/vaibhav/schulze_project/gb/IDp03_aimsgbCK/00-minimize/01-heat-723/02-diffusion",
               "/nfshome/deshmukh/vaibhav/schulze_project/gb/IDp03_aimsgbCK/00-minimize/01-heat-773/02-diffusion"]

# Chosing half-width of of grain boundary w in Angstrom
w = 10 # Check with Ovito!



for element in input_files:
    # Data import
    infile = element
    for i,j in zip(range(500,3000,500), range(0,2500,500)):    
        os.chdir = infile
        pipeline = import_file(f"{infile}/structure.dump.*")
        print(infile)

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

        def modify1(frame: int, data: DataCollection):

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

        pipeline.modifiers.append(modify1)


        # Displacement vectors

        pipeline.modifiers.append(CalculateDisplacementsModifier(reference_frame = j, minimum_image_convention = False))


        # Calculate 'Relative Displacement' of host matrix

        def modify2(frame, data, type_id_Na = 4):
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
        pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 4))

        # Calculate MSD, netMSD for x-,y-,z-direction

        #pipeline.modifiers.append(SliceModifier(
        #        distance = 0.5, 
        #        normal = (0.0, 0.0, 1.0), 
        #        slab_width = 10.0, 
        #        miller = True))

        def modify3(frame: int, data: DataCollection, w = w, type_id_Na = 4):

            #Select GB: Choice has to be double-checked in Ovito!
            #pos = data.particles.positions
            #Selection = data.particles_.create_property("Selection")
            #Selection[(pos[:,2] < data.cell[2][3] + w) ] = 1
            #Selection[(pos[:,2] > data.cell[2][3] + data.cell[2][2]/2. - w) & (pos[:,2] < data.cell[2][3] + data.cell[2][2]/2. + w) ] = 1
            #Selection[(pos[:,2] > data.cell[2][3] + data.cell[2][2] - w) ] = 1
        
            #data.apply(FreezePropertyModifier(source_property = 'Selection', 
            #                          destination_property = 'Select0',
            #                          freeze_at = 0))

            ptype = data.particles.particle_type
            ##only Na atoms
            #Selection = data.particles["Select0"]
        
            time = frame*2000*2/1000
            
            box_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]         # NEWMS

            # Ionic diffusion for box (COM corrected)
            net_x, net_y, net_z = np.sum(np.abs(box_rel), axis = 0)
            data.attributes["netMSDx_box_com"] = net_x**2 / len(box_rel) #netMSD = ionic diffusivity
            data.attributes["netMSDy_box_com"] = net_y**2 / len(box_rel)
            data.attributes["netMSDz_box_com"] = net_z**2 / len(box_rel)
            data.attributes["netMSDav_box_com"] = (net_x**2)+(net_y**2)+(net_z**2) / len(box_rel)
            
            data.attributes["Timestep"] = time
        
        pipeline.modifiers.append(modify3)

        export_file(pipeline, f"{infile}/na_diffusion_tilt_auto_new_ch_%s.txt"%(j), "txt/attr", columns=["Timestep","netMSDx_box_com", "netMSDy_box_com", "netMSDz_box_com","netMSDav_box_com"],multiple_frames=True, start_frame=j, end_frame=i)
    
