from fileinput import filename
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
import numpy as np 

def smooth_msd(path_line='dump_nvt_prod.out',file_2='msd_na.txt',step=1.0, dump_write=100, step_ovito=10,m=0,n=1000):


    pipeline = import_file(path_line)
    
    n_frames = pipeline.source.num_frames
    Na_type_id = 2

    # Simulation Parameters
    POTIM = step                  # Time step in fs
    NBLOCK = dump_write            # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
    # MSD Evaluation
    Stepsize = step_ovito 

    # set element 
    pipeline.add_to_scene()
    def setup_particle_types(frame, data):
        types = data.particles_.particle_types_
        types.type_by_id_(1).name = "Li"
        types.type_by_id_(2).name = "Ni"
        types.type_by_id_(3).name = "O"
    pipeline.modifiers.append(setup_particle_types)

    def modify1(frame: int, data: DataCollection):

        charge_dict = {"Li": 2.4, "Ni": 2.4, "O":-1.2}     # Padone charges 
        mass_dict = {"Li": 6.94,
                        "Ni": 58.6934,
                        "O": 15.9994}
        charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
        mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
        ptypes = data.particles_.particle_types_
        for i in range(data.particles.count):
            charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
            mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
    pipeline.modifiers.append(modify1)

    displmod = CalculateDisplacementsModifier(reference_frame = m, minimum_image_convention = False)#,use_frame_offset = True)
    pipeline.modifiers.append(displmod)

    def modify2(frame, data, type_id_Na = 1):
            
        #with respect to center of mass of host matrix   
        #basis = data.cell[0:3,0:3] 
        #print(basis.T)                                                                                                                                                                               
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
    #    crystal = np.dot(np.linalg.inv(basis.T).T, (displ - com_displ_host_matrix).T).T                 # Basis transformation to fractional coordinates
    #    data.particles_.create_property('Relative Displacement', data = np.dot((basis.T).T, crystal.T).T) # Basis transformation to cartesian coordinates
        data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
   # 
    import functools
    pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 1))



    def modify3(frame: int, data: DataCollection, type_id_Na = 1):

        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_x = np.mean(bulk_displ_rel[:, 0]**2, axis=0)
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_y = np.mean(bulk_displ_rel[:, 1]**2, axis=0)
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_z = np.mean(bulk_displ_rel[:, 2]**2, axis=0)
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_xy = np.mean(bulk_displ_rel[:, 0] * bulk_displ_rel[:, 1], axis=0)
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_xz = np.mean(bulk_displ_rel[:, 0] * bulk_displ_rel[:, 2], axis=0)
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_yz = np.mean(bulk_displ_rel[:, 1] * bulk_displ_rel[:, 2], axis=0)      
    
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_x = (np.sum(bulk_displ_rel[:, 0], axis = 0))**2 / len(bulk_displ_rel)
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_y = np.sum(bulk_displ_rel[:, 1], axis = 0)**2 / len(bulk_displ_rel)
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_z = np.sum(bulk_displ_rel[:, 2], axis = 0)**2 / len(bulk_displ_rel)
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_xy = np.sum(bulk_displ_rel[:, 0]) * np.sum(bulk_displ_rel[:, 1]) / len(bulk_displ_rel)
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_xz = np.sum(bulk_displ_rel[:, 0]) * np.sum(bulk_displ_rel[:, 2]) / len(bulk_displ_rel)
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_yz = np.sum(bulk_displ_rel[:, 1]) * np.sum(bulk_displ_rel[:, 2]) / len(bulk_displ_rel)
        MSDav_bulk_com = msd_x + msd_y + msd_z      
        netMSDav_bulk_com = net_x + net_y + net_z

        data.attributes['MSDx_bulk_com'] = msd_x 
        data.attributes['MSDy_bulk_com'] = msd_y
        data.attributes['MSDz_bulk_com'] = msd_z
        data.attributes['MSDxy_bulk_com'] = msd_xy
        data.attributes['MSDxz_bulk_com'] = msd_xz
        data.attributes['MSDyz_bulk_com'] = msd_yz
        data.attributes['netMSDx_bulk_com'] = net_x
        data.attributes['netMSDy_bulk_com'] = net_y
        data.attributes['netMSDz_bulk_com'] = net_z
        data.attributes['netMSDxy_bulk_com'] = net_xy
        data.attributes['netMSDxz_bulk_com'] = net_xz
        data.attributes['netMSDyz_bulk_com'] = net_yz
        data.attributes['MSDav_bulk_com'] = MSDav_bulk_com
        data.attributes['netMSDav_bulk_com'] = netMSDav_bulk_com 

        data.attributes["Timestep"] = frame*NBLOCK*POTIM/1000

        print('Timestep: ', frame*NBLOCK*POTIM/1000, 'ps')

    pipeline.modifiers.append(modify3)

    export_file(pipeline, file_2, "txt/attr",
    columns=["Timestep", "MSDx_bulk_com", "MSDy_bulk_com", "MSDz_bulk_com", "MSDxy_bulk_com", "MSDxz_bulk_com", "MSDyz_bulk_com", "netMSDx_bulk_com", "netMSDy_bulk_com", "netMSDz_bulk_com", "netMSDxy_bulk_com", "netMSDxz_bulk_com", "netMSDyz_bulk_com", "MSDav_bulk_com", "netMSDav_bulk_com"],multiple_frames=True, start_frame=m, end_frame=n)


print('Computation for MSD has started, please wait...')

import tarfile
import os

li =  [800]
no = [4,4,4,4,4,4]

#file_paths = ["/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_523k_4_hdf5/hena_1_struct_eq_big_523k_4/hena_1_struct_eq_big_523k_4.tar.bz2",
#               "/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_573k_4_hdf5/hena_1_struct_eq_big_573k_4/hena_1_struct_eq_big_573k_4.tar.bz2",
#               "/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_623k_4_hdf5/hena_1_struct_eq_big_623k_4/hena_1_struct_eq_big_623k_4.tar.bz2",
#               "/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_673k_4_hdf5/hena_1_struct_eq_big_673k_4/hena_1_struct_eq_big_673k_4.tar.bz2",
#               "/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_723k_4_hdf5/hena_1_struct_eq_big_723k_4/hena_1_struct_eq_big_723k_4.tar.bz2",
#               "/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_773k_4_hdf5/hena_1_struct_eq_big_773k_4/hena_1_struct_eq_big_773k_4.tar.bz2",]

#import os
#import tarfile
#import bz2

#for file_path in file_paths:
#    folder_path = os.path.dirname(file_path)
#    output_file_name = 'dump_nvt_prod.out'#

#    try:
#        if file_path.endswith(".bz2"):
#            with open(output_file_name, 'wb') as output_file, bz2.BZ2File(file_path, 'rb') as input_file:
#                output_file.write(input_file.read())
#        elif file_path.endswith(".bz"):
#            with tarfile.open(file_path, 'r:bz2') as tar:
#                tar.extract(output_file_name, path=folder_path)
#        print(f"File extraction completed for {file_path}.")
#    except (tarfile.TarError, IOError) as e:
#        print(f"An error occurred while extracting the file {file_path}: {str(e)}")
#    # Additional processing or actions with the extracted files can be done here

input_files = ["/nfshome/sadowski/work/LiNiO2_data_base_Sabrina/08_MD_for_AL/05_long_MD_test_for_statistics/800K/dump_nvt_prod.out"]

for i,j,k in zip(li,input_files,no):
    pipeline = import_file(j)
    frame = pipeline.source.num_frames
    print(frame)
    d = 100 
    for e in range(0,pipeline.source.num_frames,pipeline.source.num_frames//d):
        #if e+pipeline.source.num_frames//100 > pipeline.source.num_frames:
        #    break
        if e == 0:
            m = 0
            n = pipeline.source.num_frames//d 
        
        else:
        #    #print(e, e+pipeline.source.num_frames//d)
            m = e-pipeline.source.num_frames//d
            n = e 
            
            print(m,n)
            #for m,n in zip(range(500,3000,500), range(0,2500,500)):
            print('------------------------------------------------------------------------------------------------------------------')
            print('------------------------------------------MSD_Calculations Info---------------------------------------------------')
            print('------------------------------------[Files (Temperature/Location/No)]---------------------------------------------')
            print(i,j,k)
            print('-------------------------------------------Types of Calulations---------------------------------------------------')
            print('1) MSD of Na ions, x, y, z, xy, xz, yz, av, net_av')
            print('2) Chharge MSD of Na ions, x, y, z, xy, xz, yz, av, net_av')
            print('------------------------------------------------------------------------------------------------------------------')
            smooth_msd(path_line=j,file_2='msd_%s_%s_%s_%s.txt'%(i,k,m,n),step=1.00, dump_write=5000, step_ovito=1, m=m, n=n)
            print('-----------------------------------------------Next calc----------------------------------------------------------')
            print('---------------------------------------------------------------------------------------------------------------------')
            print('Done with MSD calculations, please check the files for the results. Thank you for your patience.')
