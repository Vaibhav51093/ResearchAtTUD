from fileinput import filename
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
import numpy as np 

def diff_4(path_line='dump_nvt_prod.out',file_2='msd_na.txt',step=1.0, dump_write=100, step_ovito=10):

    file = path_line

    # Data import 
    pipeline = import_file(file)

    # Simulation Parameters
    POTIM = step                  # Time step in fs
    NBLOCK = dump_write            # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
    # MSD Evaluation
    Stepsize = step_ovito 

    # Print the list of input particle type
    pipeline.modifiers.append(SmoothTrajectoryModifier(minimum_image_convention=True,window_size=Stepsize))

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


    # Displacement vectors

    pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

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

    def modify3(frame: int, data: DataCollection, type_id_Na = 2):

        ptype = data.particles.particle_type

        time = frame*NBLOCK*POTIM/1000

        pos = data.particles['Position'][(ptype == type_id_Na)]
        #print(len(pos[:,0]))
        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        bulk_displ = data.particles['Displacement'][(ptype == type_id_Na)]
        #print(bulk_displ_rel)
        # Note:- MSD are the absolute displacement suqared, therefore we need absolute displacement 
        #        The reason to have negative MSD in xy plane  
        #msd_x, msd_y, msd_z = np.mean(np.abs(bulk_displ_rel)**2, axis = 0)
        #data.attributes["MSDx_bulk_com"] = msd_x  # com = center of mass /displacements are corrected for com shift of host matrix
        #data.attributes["MSDy_bulk_com"] = msd_y  # MSD = tracer diffucivity
        #data.attributes["MSDz_bulk_com"] = msd_z
        #data.attributes["MSDav_bulk_com"] = msd_x+msd_y+msd_z
        
        # The follwoing code implement cros_displacement algo to get offdiagonal compoenets of MSD 
        # Basically I am caulcating the varienece 
        #sum_of_pos_xy = np.mean((np.abs(bulk_displ_rel[:,0]) + np.abs(bulk_displ_rel[:,1]))**2)
        #sum_of_pos_xz = np.mean((np.abs(bulk_displ_rel[:,0]) + np.abs(bulk_displ_rel[:,2]))**2)
        #sum_of_pos_yz = np.mean((np.abs(bulk_displ_rel[:,1]) + np.abs(bulk_displ_rel[:,2]))**2)
        
        #cd_xy = cd[0][1]
        #cd_xz = cd[0][2]
        #cd_yz = cd[1][2]
        
        #print(cd_xy)
        
        #data.attributes["MSDxy_bulk_com"] = 0.5*(sum_of_pos_xy - msd_x - msd_y) #np.mean(cd_xy**2) #np.mean(np.array(cd_xy)**2, axis=1) 
        #data.attributes["MSDxz_bulk_com"] = 0.5*(sum_of_pos_xz - msd_x - msd_z) #np.mean(cd_xz**2) #np.mean(np.array(cd_xz)**2, axis=1) 
        #data.attributes["MSDyz_bulk_com"] = 0.5*(sum_of_pos_yz - msd_y - msd_z) #np.mean(cd_yz**2) #np.mean(np.array(cd_yz)**2, axis=1)

        net_x, net_y, net_z = np.sum(np.abs(bulk_displ_rel), axis = 0)
        data.attributes["netMSDx_bulk_com"] = net_x**2 / len(bulk_displ_rel) #netMSD = ionic diffusivity
        data.attributes["netMSDy_bulk_com"] = net_y**2 / len(bulk_displ_rel)
        data.attributes["netMSDz_bulk_com"] = net_z**2 / len(bulk_displ_rel)
        data.attributes["netMSDav_bulk_com"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ_rel)
        
        # The follwoing code implement cros_displacement algo to get offdiagonal compoenets of MSD 
        #sum_of_pos_xy_n = ((np.sum(np.abs(bulk_displ_rel[:,0])) + np.abs(np.sum(bulk_displ_rel[:,1])))**2) / len(bulk_displ_rel)
        #sum_of_pos_xz_n = ((np.sum(np.abs(bulk_displ_rel[:,0])) + np.abs(np.sum(bulk_displ_rel[:,2])))**2) / len(bulk_displ_rel)
        #sum_of_pos_yz_n = ((np.sum(np.abs(bulk_displ_rel[:,1])) + np.abs(np.sum(bulk_displ_rel[:,2])))**2) / len(bulk_displ_rel)
        
        #data.attributes["netMSDxy_bulk_com"] = 0.5*(sum_of_pos_xy_n - net_x**2/len(np.abs(bulk_displ_rel)) - net_y**2/len(np.abs(bulk_displ_rel))) #np.mean(cd_xy**2) #np.mean(np.array(cd_xy)**2, axis=1) 
        #data.attributes["netMSDxz_bulk_com"] = 0.5*(sum_of_pos_xz_n - net_x**2/len(np.abs(bulk_displ_rel)) - net_z**2/len(np.abs(bulk_displ_rel))) #np.mean(cd_xz**2) #np.mean(np.array(cd_xz)**2, axis=1) 
        #data.attributes["netMSDyz_bulk_com"] = 0.5*(sum_of_pos_yz_n - net_y**2/len(np.abs(bulk_displ_rel)) - net_z**2/len(np.abs(bulk_displ_rel)))
        
        data.attributes["Timestep"] = time
     
    pipeline.modifiers.append(modify3)

    export_file(pipeline, file_2, "txt/attr",
    columns=["Timestep", "netMSDx_bulk_com", "netMSDy_bulk_com", "netMSDz_bulk_com","netMSDav_bulk_com"],
    multiple_frames=True)
    
    

