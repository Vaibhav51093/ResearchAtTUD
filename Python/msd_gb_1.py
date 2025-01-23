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
    #pipeline.modifiers.append(WrapPeriodicImagesModifier())
    
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

    displmod = CalculateDisplacementsModifier(reference_frame=m, minimum_image_convention=False)
    pipeline.modifiers.append(displmod)

    # Define region using SliceModifier
    region_center = 0.5  # Center of the region along y-axis (assuming normalized coordinates)
    region_width = 20.0  # Width of the region to select

    dist = shift*3
    gb_width = 10.0
    
    # Slice:
    #pipeline.modifiers.append(SliceModifier(
    #    distance = dist, 
    #    normal = (0.0, 1.0, 0.0), 
    #    slab_width = gb_width, 
    #    inverse = True, 
    #    select = True))
    
    pipeline.modifiers.append(SliceModifier(
        distance = 0.5,
        normal = (0.0, 1.0, 0.0),
        slab_width = gb_width,
        inverse = True,
        select = True,
        miller = True))
    
    # Invert selection:
    pipeline.modifiers.append(InvertSelectionModifier())

    # Delete selected:
    pipeline.modifiers.append(DeleteSelectedModifier())
    
    # Select type:
    pipeline.modifiers.append(SelectTypeModifier(types = {2}))
    
    # unwrap the trajectory 
    # Unwrap trajectories:
    #pipeline.modifiers.append(UnwrapTrajectoriesModifier())
    
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
    
    # Invert selection:
    #pipeline.modifiers.append(InvertSelectionModifier())

    # Delete selected:
    #pipeline.modifiers.append(DeleteSelectedModifier())

    def modify3(frame, data: DataCollection, type_id_Na=2):
        ptype = data.particles['Particle Type']
        selected = data.particles['Selection'] > 0
        displ_rel = data.particles['Relative Displacement'][np.logical_and(ptype == type_id_Na, selected)]

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

        print('Timestep: ', frame * NBLOCK * POTIM / 1000, 'ps')

    pipeline.modifiers.append(modify3)

    export_file(pipeline, file_2, "txt/attr",
    columns=["Timestep", "MSDx", "MSDy", "MSDz", "MSDxy", "MSDxz", "MSDyz", "netMSDx", "netMSDy", "netMSDz", "netMSDxy", "netMSDxz", "netMSDyz", "MSDav", "netMSDav"], multiple_frames=True, start_frame=m, end_frame=n)

print('Computation for MSD has started, please wait...')



import tarfile
import os

li = [523,573,623,673,723,773]
no = [6,6,6,6,6,6]

input_files = ["/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_13_120/01-heat-523/02-diffusion/structure.dump.*.gz",
               "/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_13_120/01-heat-573/02-diffusion/structure.dump.*.gz",
               "/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_13_120/01-heat-623/02-diffusion/structure.dump.*.gz",
               "/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_13_120/01-heat-673/02-diffusion/structure.dump.*.gz",
               "/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_13_120/01-heat-723/02-diffusion/structure.dump.*.gz",
               "/work/scratch/vd80naku/pyiron/projects/remote/NASICON/gb_new/gb_13_120/01-heat-773/02-diffusion/structure.dump.*.gz"]
       
for i,j,k in zip(li,input_files,no):
    pipeline = import_file(j)
    frame = pipeline.source.num_frames
    d = 10 
    
    #common_na = compute_common_ids(path_line=j,step=2.0, dump_write=200, step_ovito=1, m=0, n=pipeline.source.num_frame

    for f in range(0, pipeline.source.num_frames, 100):

        for v in range(0+f, pipeline.source.num_frames, 450):

            m = v
            n = v+500  # 2 ns interval 
            if n > pipeline.source.num_frames:
                break
            else:
                #for m,n in zip(range(500,3000,500), range(0,2500,500)):
                print('------------------------------------------------------------------------------------------------------------------')
                print('------------------------------------------MSD_Calculations Info---------------------------------------------------')
                print('------------------------------------[Files (Temperature/Location/No)]---------------------------------------------')
                print(i,j,k)
                print('-------------------------------------------Types of Calulations---------------------------------------------------')
                print('1) MSD of Na ions, x, y, z, xy, xz, yz, av, net_av')
                print('2) Chharge MSD of Na ions, x, y, z, xy, xz, yz, av, net_av')
                print('------------------------------------------------------------------------------------------------------------------')
                #common_na = compute_common_ids(path_line=j,step=2.0, dump_write=200, step_ovito=1, m=m, n=n)
                #print(common_na)
                smooth_msd(path_line=j,file_2='msd_gb_1_%s_%s_%s_%s.txt'%(i,k,m,n),step=2.0, dump_write=200, step_ovito=1, m=m, n=n) #, common_indices=common_na)
                print('-----------------------------------------------Next calc----------------------------------------------------------')
                print('---------------------------------------------------------------------------------------------------------------------')
                print('Done with MSD calculations, please check the files for the results. Thank you for your patience.')

