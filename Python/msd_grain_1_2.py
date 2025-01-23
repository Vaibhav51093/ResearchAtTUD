from fileinput import filename
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
import numpy as np 
import functools

def smooth_msd(path_line='dump_nvt_prod.out', file_2='msd_na.txt', step=1.0, dump_write=100, step_ovito=10, m=0, n=1000):

    pipeline = import_file(path_line)
    
    n_frames = pipeline.source.num_frames
    Na_type_id = 2

    # Simulation Parameters
    POTIM = step                  # Time step in fs
    NBLOCK = dump_write           # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
    # MSD Evaluation
    Stepsize = step_ovito 

    def setup_particle_types(frame, data):
        types = data.particles_.particle_types_
        types.type_by_id_(1).name = "O"
        types.type_by_id_(2).name = "Na"
        types.type_by_id_(3).name = "Zr"
        types.type_by_id_(4).name = "Si"
        types.type_by_id_(5).name = "P"
        
    pipeline.modifiers.append(setup_particle_types)
    
    data_initial = pipeline.compute(frame=0)
    
    Z_Length = data_initial.cell[1][1]
    
    shift = Z_Length/4
        
    # Affine transformation:
    #pipeline.modifiers.append(AffineTransformationModifier(
    #    transformation = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, shift], [0.0, 0.0, 1.0, 0.0]], 
    #    operate_on = {'dislocations', 'surfaces', 'particles', 'voxels', 'vector_properties'}))

    # Wrap at periodic boundaries:
    pipeline.modifiers.append(WrapPeriodicImagesModifier())
    
    # Set charges and mass
    def modify1(frame, data: DataCollection):
        charge_dict = { "Si": 2.4, "Zr": 2.4, "P": 3.0, "O": -1.2, "Na": 0.6 }
        mass_dict = {"Zr": 91.224, "Si": 28.086, "P": 30.974, "O": 15.999, "Na": 22.99 }
        charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
        mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
        ptypes = data.particles_.particle_types_
        for i in range(data.particles.count):
            charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
            mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]

    pipeline.modifiers.append(modify1)

    displmod = CalculateDisplacementsModifier(reference_frame=m, minimum_image_convention=True)
    pipeline.modifiers.append(displmod)

    # Define region using SliceModifier
    #region_center = 0.5  # Center of the region along y-axis (assuming normalized coordinates)
    #region_width = 20.0  # Width of the region to select

    dist = shift*3
    gb_width = (Z_Length/2) - 10.0
    
    # Slice:
    #pipeline.modifiers.append(SliceModifier(
    #    distance = dist, 
    #    normal = (0.0, 1.0, 0.0), 
    #    slab_width = gb_width, 
    #    inverse = True, 
    #    select = True))
    
    pipeline.modifiers.append(SliceModifier(
        distance = 0.75,
        normal = (0.0, 1.0, 0.0),
        slab_width = gb_width,
        inverse = True,
        select = True,
        miller = True))

    def modify2(frame, data, type_id_Na=2):
        ptypes = data.particles['Particle Type']
        displ_matrix = data.particles['Displacement'][(ptypes != type_id_Na)]
        displ = data.particles['Displacement']
        masses = data.particles['Mass'][(ptypes != type_id_Na)]

        # Calculate center of mass displacement of host matrix
        total_mass_matrix = np.sum(masses)
        com_displ_host_matrix = np.sum(displ_matrix * masses[:, None], axis=0) / total_mass_matrix
        data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
        data.particles_.create_property('Relative Displacement', data=displ - com_displ_host_matrix)

    pipeline.modifiers.append(functools.partial(modify2, type_id_Na=2))
    
    # Compute particle in each frame and create a list for each frame and get common index for all frames
    from functools import reduce
    
    na_indices = []
    
    for i in range(0, n_frames, Stepsize):
        data_0 = pipeline.compute(frame=i)
        selected = data_0.particles['Selection'] > 0
        index = data_0.particles['Particle Identifier'][np.logical_and(data_0.particles['Particle Type'] == Na_type_id, selected)]
        na_indices.append(index)
        #print(index)
        #print('Frame: ', i, 'Na ions: ', len(index))
        
    # Convert lists to sets
    sets = [set(lst) for lst in na_indices]
    
    # Find common indices using reduce and set intersection
    common_indices = list(reduce(lambda x, y: x.intersection(y), sets))
    print(common_indices)
    print('Common Na ions: ', len(common_indices))

    def modify3(frame, data: DataCollection, type_id_Na=2):
        # old
        ptype = data.particles['Particle Type']
        selected = data.particles['Selection'] > 1
        displ_rel = data.particles['Relative Displacement'][(common_indices )]
        
              
        msd_x = np.mean(displ_rel[:, 0]**2, axis=0)
        msd_y = np.mean(displ_rel[:, 1]**2, axis=0)
        msd_z = np.mean(displ_rel[:, 2]**2, axis=0)
        msd_xy = np.mean(displ_rel[:, 0] * displ_rel[:, 1], axis=0)
        msd_xz = np.mean(displ_rel[:, 0] * displ_rel[:, 2], axis=0)
        msd_yz = np.mean(displ_rel[:, 1] * displ_rel[:, 2], axis=0)

        net_x = (np.sum(displ_rel[:, 0], axis=0))**2 / len(displ_rel)
        net_y = (np.sum(displ_rel[:, 1], axis=0))**2 / len(displ_rel)
        net_z = (np.sum(displ_rel[:, 2], axis=0))**2 / len(displ_rel)
        net_xy = np.sum(displ_rel[:, 0]) * np.sum(displ_rel[:, 1]) / len(displ_rel)
        net_xz = np.sum(displ_rel[:, 0]) * np.sum(displ_rel[:, 2]) / len(displ_rel)
        net_yz = np.sum(displ_rel[:, 1]) * np.sum(displ_rel[:, 2]) / len(displ_rel)

        MSDav = msd_x + msd_y + msd_z
        netMSDav = net_x + net_y + net_z

        data.attributes['MSDx'] = msd_x 
        data.attributes['MSDy'] = msd_y
        data.attributes['MSDz'] = msd_z
        data.attributes['MSDxy'] = msd_xy
        data.attributes['MSDxz'] = msd_xz
        data.attributes['MSDyz'] = msd_yz
        data.attributes['netMSDx'] = net_x
        data.attributes['netMSDy'] = net_y
        data.attributes['netMSDz'] = net_z
        data.attributes['netMSDxy'] = net_xy
        data.attributes['netMSDxz'] = net_xz
        data.attributes['netMSDyz'] = net_yz
        data.attributes['MSDav'] = MSDav
        data.attributes['netMSDav'] = netMSDav 

        data.attributes["Timestep"] = frame * NBLOCK * POTIM / 1000

        #print('Timestep: ', frame * NBLOCK * POTIM / 1000, 'ps')

    pipeline.modifiers.append(modify3)

    export_file(pipeline, file_2, "txt/attr",
    columns=["Timestep", "MSDx", "MSDy", "MSDz", "MSDxy", "MSDxz", "MSDyz", "netMSDx", "netMSDy", "netMSDz", "netMSDxy", "netMSDxz", "netMSDyz", "MSDav", "netMSDav"], multiple_frames=True, start_frame=m, end_frame=n)

print('Computation for MSD has started, please wait...')



import tarfile
import os

li =  [523,573,623,673,723,773]
no = [6,6,6,6,6,6]

#file_paths = ["/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_523k_43_hdf5/hena_1_struct_eq_big_523k_43/dump_nvt_prod.out.tar.gz",
#               "/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_573k_43_hdf5/hena_1_struct_eq_big_573k_43/dump_nvt_prod.out.tar.gz",
#               "/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_623k_43_hdf5/hena_1_struct_eq_big_623k_43/dump_nvt_prod.out.tar.gz",
#               "/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_673k_43_hdf5/hena_1_struct_eq_big_673k_43/dump_nvt_prod.out.tar.gz",
#               "/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_723k_43_hdf5/hena_1_struct_eq_big_723k_43/dump_nvt_prod.out.tar.gz",
#               "/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_773k_43_hdf5/hena_1_struct_eq_big_773k_43/dump_nvt_prod.out.tar.gz",]

#for file_path in file_paths:
#    folder_path = os.path.dirname(file_path)
#    output_file_name = 'dump_nvt_prod.out'##

#    try:
#        with tarfile.open(file_path, 'r:gz') as tar:
#            tar.extract(output_file_name, path=folder_path)
#        print(f"File extraction completed for {file_path}.")
#    except tarfile.TarError as e:
#        print(f"An error occurred while extracting the file {file_path}:", str(e))
#gb_13_230/relaxation_2/gb_13_120
input_files = ["/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_5_120/01-heat-523/02-diffusion/structure.dump.*.gz",
               "/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_5_120/01-heat-573/02-diffusion/structure.dump.*.gz",
               "/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_5_120/01-heat-623/02-diffusion/structure.dump.*.gz",
               "/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_5_120/01-heat-673/02-diffusion/structure.dump.*.gz",
               "/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_5_120/01-heat-723/02-diffusion/structure.dump.*.gz",
               "/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_5_120/01-heat-773/02-diffusion/structure.dump.*.gz"]
       
for i,j,k in zip(li,input_files,no):
    pipeline = import_file(j)
    frame = pipeline.source.num_frames
    d = 10 
    for e in range(0,pipeline.source.num_frames,pipeline.source.num_frames//d):
        if e+pipeline.source.num_frames//5 > pipeline.source.num_frames:
            break
        elif pipeline.source.num_frames == pipeline.source.num_frames//d:
            m = 0
            n = pipeline.source.num_frames
        else:
            #print(e, e+pipeline.source.num_frames//d)
            m = e
            n = e+pipeline.source.num_frames//5
            
            #for m,n in zip(range(500,3000,500), range(0,2500,500)):
            print('------------------------------------------------------------------------------------------------------------------')
            print('------------------------------------------MSD_Calculations Info---------------------------------------------------')
            print('------------------------------------[Files (Temperature/Location/No)]---------------------------------------------')
            print(i,j,k)
            print('-------------------------------------------Types of Calulations---------------------------------------------------')
            print('1) MSD of Na ions, x, y, z, xy, xz, yz, av, net_av')
            print('2) Chharge MSD of Na ions, x, y, z, xy, xz, yz, av, net_av')
            print('------------------------------------------------------------------------------------------------------------------')
            smooth_msd(path_line=j,file_2='msd_grain_1_%s_%s_%s_%s.txt'%(i,k,m,n),step=2.0, dump_write=200, step_ovito=1, m=m, n=n)
            print('-----------------------------------------------Next calc----------------------------------------------------------')
            print('---------------------------------------------------------------------------------------------------------------------')
            print('Done with MSD calculations, please check the files for the results. Thank you for your patience.')

