from fileinput import filename
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
import numpy as np 
import tidynamics

def smooth_msd(path_line='dump_nvt_prod.out',file_2='msd_na.txt',step=1.0, dump_write=100, step_ovito=10):


    pipeline = import_file(path_line)
    n_frames = pipeline.source.num_frames
    Na_type_id = 2

    dt = step * dump_write 
    half_frame = pipeline.source.num_frames #// 2

    every_nth_frame = step_ovito

    # set element 
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



    displmod = CalculateDisplacementsModifier(minimum_image_convention = False,use_frame_offset = True)
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


    def MSDx_bulk_com(data, type_id_Na = 2):
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_x = np.mean(bulk_displ_rel[:, 0]**2, axis=0)
        return msd_x
    
    def MSDy_bulk_com(data, type_id_Na = 2):
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_y = np.mean(bulk_displ_rel[:, 1]**2, axis=0)
        return msd_y
    
    def MSDz_bulk_com(data, type_id_Na = 2):
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_z = np.mean(bulk_displ_rel[:, 2]**2, axis=0)
        return msd_z
    
    def MSDxy_bulk_com(data, type_id_Na = 2):
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_xy = np.mean(bulk_displ_rel[:, 0] * bulk_displ_rel[:, 1], axis=0)
        return msd_xy
    
    def MSDxz_bulk_com(data, type_id_Na = 2):       
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_xz = np.mean(bulk_displ_rel[:, 0] * bulk_displ_rel[:, 2], axis=0)
        return msd_xz
    
    def MSDyz_bulk_com(data, type_id_Na = 2):
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        msd_yz = np.mean(bulk_displ_rel[:, 1] * bulk_displ_rel[:, 2], axis=0)
        return msd_yz       
    
    def netMSDx_bulk_com(data, type_id_Na = 2):
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_x = np.sum(bulk_displ_rel[:, 0], axis = 0)
        return net_x**2 / len(bulk_displ_rel)
    
    def netMSDy_bulk_com(data, type_id_Na = 2):     
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_y = np.sum(bulk_displ_rel[:, 1], axis = 0)
        return net_y**2 / len(bulk_displ_rel)
    
    def netMSDz_bulk_com(data, type_id_Na = 2):
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_z = np.sum(bulk_displ_rel[:, 2], axis = 0)
        return net_z**2 / len(bulk_displ_rel)
    
    def netMSDxy_bulk_com(data, type_id_Na = 2):
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_xy = np.sum(bulk_displ_rel[:, 0]) * np.sum(bulk_displ_rel[:, 1]) / len(bulk_displ_rel)
        return net_xy
    
    def netMSDxz_bulk_com(data, type_id_Na = 2):
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_xz = np.sum(bulk_displ_rel[:, 0]) * np.sum(bulk_displ_rel[:, 2]) / len(bulk_displ_rel)
        return net_xz   
    
    def netMSDyz_bulk_com(data, type_id_Na = 2):        
        ptype = data.particles.particle_type
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]        
        net_yz = np.sum(bulk_displ_rel[:, 1]) * np.sum(bulk_displ_rel[:, 2]) / len(bulk_displ_rel)
        return net_yz
    
    def MSDav_bulk_com(data, type_id_Na = 2):
        return MSDx_bulk_com(data, type_id_Na) + MSDy_bulk_com(data, type_id_Na) + MSDz_bulk_com(data, type_id_Na)      
    
    def netMSDav_bulk_com(data, type_id_Na = 2):
        return netMSDx_bulk_com(data, type_id_Na) + netMSDy_bulk_com(data, type_id_Na) + netMSDz_bulk_com(data, type_id_Na)
    
    # Smooth the MSD values

    #msd_bulk = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    #msd_bulk_net = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

    #for i in range(1, half_frame, every_nth_frame):
    #    print(f'frame = {i}')
    #    displ.frame_offset = -i 
    #    for frame in range(i, n_frames):
    #        data = pipeline.compute(frame)
    #         msd_bulk += np.array((MSDx_bulk_com(data, Na_type_id), MSDy_bulk_com(data, Na_type_id), MSDz_bulk_com(data, Na_type_id), MSDav_bulk_com(data, Na_type_id), MSDxy_bulk_com(data, Na_type_id), MSDxz_bulk_com(data, Na_type_id), MSDyz_bulk_com(data, Na_type_id)))
    #         msd_bulk_net += np.array((netMSDx_bulk_com(data, Na_type_id), netMSDy_bulk_com(data, Na_type_id), netMSDz_bulk_com(data, Na_type_id), netMSDav_bulk_com(data, Na_type_id), netMSDxy_bulk_com(data, Na_type_id), netMSDxz_bulk_com(data, Na_type_id), netMSDyz_bulk_com(data, Na_type_id)))
    #    msd_bulk = np.concatenate((msd_bulk, (msd_bulk/(n_frames - i))[np.newaxis, :]), axis=0)
    #    msd_bulk_net = np.concatenate((msd_bulk_net, (msd_bulk_net/(n_frames - i))[np.newaxis, :]), axis=0)

    #msd_bulk = np.empty((0, 7))
    #msd_bulk_net = np.empty((0, 7))



    #for i in range(1, half_frame, every_nth_frame):
    #    print(f'frame = {i}')
    #    
    #    #displ.frame_offset = -i#

    #    msd_bulk_temp = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    #    msd_bulk_net_temp = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    #    for frame in range(i, n_frames):
    #        data = pipeline.compute(frame)
    #        msd_bulk_temp += np.array((MSDx_bulk_com(data, Na_type_id), MSDy_bulk_com(data, Na_type_id), MSDz_bulk_com(data, Na_type_id), MSDav_bulk_com(data, Na_type_id), MSDxy_bulk_com(data, Na_type_id), MSDxz_bulk_com(data, Na_type_id), MSDyz_bulk_com(data, Na_type_id)))
    #        msd_bulk_net_temp += np.array((netMSDx_bulk_com(data, Na_type_id), netMSDy_bulk_com(data, Na_type_id), netMSDz_bulk_com(data, Na_type_id), netMSDav_bulk_com(data, Na_type_id), netMSDxy_bulk_com(data, Na_type_id), netMSDxz_bulk_com(data, Na_type_id), netMSDyz_bulk_com(data, Na_type_id)))#

    #    msd_bulk = np.concatenate((msd_bulk, (msd_bulk_temp / (n_frames - i))[np.newaxis, :]), axis=0)
    #    msd_bulk_net = np.concatenate((msd_bulk_net, (msd_bulk_net_temp / (n_frames - i))[np.newaxis, :]), axis=0)

        
    # Write the smoothed MSD values to the output file
    #t = np.arange(len(msd_bulk), dtype=float)
    #t*=(dt*every_nth_frame)
    #np.savetxt(file_2, np.column_stack((t, msd_bulk, msd_bulk_net)), delimiter=" ", header="Timestep(ps), MSDx_bulk_com, MSDy_bulk_com, MSDz_bulk_com, MSDav_bulk_com, MSDxy_bulk_com, MSDxz_bulk_com, MSDyz_bulk_com, netMSDx_bulk_com, netMSDy_bulk_com, netMSDz_bulk_com, netMSDav_bulk_com, netMSDxy_bulk_com, netMSDxz_bulk_com, netMSDyz_bulk_com")   
    

smooth_msd(path_line="/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass/boilot_623K_2_05_cryt_523_k_ps_new_hdf5/boilot_623K_2_05_cryt_523_k_ps_new/dump_nvt_prod.out",file_2='msd_na_523_tensor_crystal_msd_smooth_2.txt',step=1.0, dump_write=100, step_ovito=1000)