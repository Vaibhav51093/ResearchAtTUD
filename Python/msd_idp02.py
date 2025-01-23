from fileinput import filename
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
import numpy as np 
#import tidynamics

def smooth_msd(path_line='dump_nvt_prod.out',file_2='msd_na.txt',step=1.0, dump_write=100, step_ovito=10,m=0,n=1000):

    pipeline = import_file(path_line)
    
    n_frames = pipeline.source.num_frames
    Na_type_id = 4

    # Simulation Parameters
    POTIM = step                  # Time step in fs
    NBLOCK = dump_write            # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
    # MSD Evaluation
    Stepsize = step_ovito 

    def setup_particle_types(frame, data):
        types = data.particles_.particle_types_
        types.type_by_id_(1).name = "O"
        types.type_by_id_(2).name = "Na"
        types.type_by_id_(3).name = "Zr"
        types.type_by_id_(4).name = "Si"
        types.type_by_id_(5).name = "P"
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

    displmod = CalculateDisplacementsModifier(reference_frame = m, minimum_image_convention = False)#,use_frame_offset = True)
    pipeline.modifiers.append(displmod)

    def modify2(frame, data, type_id_Na = 2):
            
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
    pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))



    def modify3(frame: int, data: DataCollection, type_id_Na = 2):

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

li =  [523,573,623,673,723,773]
no = [2,2,2,2,2,2]

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

input_files = ["/nfshome/deshmukh/vaibhav/schulze_project/gb/gb_new/gb_2/01_heat_523/02_diffusion/structure.dump.*",
               "/nfshome/deshmukh/vaibhav/schulze_project/gb/gb_new/gb_2/02_heat_573/02_diffusion/structure.dump.*",
               "/nfshome/deshmukh/vaibhav/schulze_project/gb/gb_new/gb_2/03_heat_623/02_diffusion/structure.dump.*",
               "/nfshome/deshmukh/vaibhav/schulze_project/gb/gb_new/gb_2/04_heat_673/02_diffusion/structure.dump.*",
               "/nfshome/deshmukh/vaibhav/schulze_project/gb/gb_new/gb_2/05_heat_723/02_diffusion/structure.dump.*",
               "/nfshome/deshmukh/vaibhav/schulze_project/gb/gb_new/gb_2/06_heat_773/02_diffusion/structure.dump.*",]


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
            smooth_msd(path_line=j,file_2='msd_na_idp_%s_%s_%s_%s.txt'%(i,k,m,n),step=2.0, dump_write=200, step_ovito=1, m=m, n=n)
            print('-----------------------------------------------Next calc----------------------------------------------------------')
            print('---------------------------------------------------------------------------------------------------------------------')
            print('Done with MSD calculations, please check the files for the results. Thank you for your patience.')



