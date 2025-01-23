    def smooth_msd_ovito(job,filename='dump_nvt_prod.out',file_2='msd_na.txt',step=0.001,dump_write=100,step_ovito=1):
        
        '''
            Smooth MSD calculations using OVITO, and has the following modules 
            1. Smooth MSD tracer with and without COM corrected 
            2. Smooth MSD charge with and without COM corrected
        '''
        # Filename and path for dump file 
        # i = reference_frame, ovito stepsize 

        import glob

        file = job.working_directory + '/%s'%(filename)
        
        pipeline = import_file(file)
        
        n_frames = pipeline.source.num_frames       # Total number of frames 
        Na_type_id = 2                              # Sodium ion type

        #tot_list = glob.glob(file)
        #list = [file for file in tot_list if int(file.split(".dump")[-1])%every_animation_frame == 0]

        # Simulation Parameters
        POTIM = step #= 0.001                 # Time step in fs
        NBLOCK = dump_write #= 100            # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
        # MSD Evaluation
        #timestep between frames in ps
        dt = POTIM*1000 # 1 fs timestep ---> 100 frames 
        #print(t, dt) 

        #Definitions for moving average
        half_frame = pipeline.source.num_frames//2                
        every_nth_frame = step_ovito
        print(half_frame)

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

            #charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
            charge_dict = {"O": -1.2, "Na": 0.6, "Zr": 2.4, "Si":2.4, "P":3.0}     # Padone charges
            mass_dict = {"O": 15.999,
                            "Na": 22.99,
                            "Zr": 91.224,
                            "Si": 28.086,
                            "P": 30.974}
            #mass_dict = {"Zr": 91.224,
            #                "Si": 28.086,
            #                "P": 30.974,
            #                "O": 15.999,
            #                "Na": 22.99}
            charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
            mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
            ptypes = data.particles_.particle_types_
            for i in range(data.particles.count):
                charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
        pipeline.modifiers.append(modify1)

        # Displacement vectors

        # Displacement vectors
        displ = CalculateDisplacementsModifier(minimum_image_convention = False, use_frame_offset = True )
        pipeline.modifiers.append(displ)

        # Calculate 'Relative Displacement' of host matrix

        def modify2(frame, data, type_id_Na = 2):
            
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
        pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

        def Calculate_MSD_grain_raw(data, type_id_Na):
            ptypes = data.particles["Particle Type"]
            #Selection = data.particles["Select0"] #####neuMS
            grain_displ = data.particles['Displacement'][(ptypes == type_id_Na)]
            return np.mean(grain_displ**2, axis=0)

        def Calculate_netMSD_grain_raw(data, type_id_Na):
            ptypes = data.particles["Particle Type"]
            #Selection = data.particles["Select0"] #####neuMS
            grain_displ = data.particles['Displacement'][(ptypes == type_id_Na)]
            return np.mean(grain_displ, axis = 0)**2

        def Calculate_MSD_grain_com(data, type_id_Na):
            ptypes = data.particles["Particle Type"]
            #Selection = data.particles["Select0"] #####neuMS
            grain_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na)]
            return np.mean(grain_displ_rel**2, axis = 0)

        def Calculate_netMSD_grain_com(data, type_id_Na):
            ptypes = data.particles["Particle Type"]
            #Selection = data.particles["Select0"] #####neuMS
            grain_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na)]
            return np.mean(grain_displ_rel, axis = 0)**2


        Na_MSD_bulk = [np.array((0.0, 0.0, 0.0))]
        net_Na_MSD_bulk = [np.array((0.0, 0.0, 0.0))]
        Na_MSD_com_bulk = [np.array((0.0, 0.0, 0.0))]
        net_Na_MSD_com_bulk = [np.array((0.0, 0.0, 0.0))]

        for i in (every_nth_frame, half_frame, every_nth_frame):
            
            print(i)

            displ.frame_offset = -i
            Na_msd_bulk = np.array((0.0, 0.0, 0.0))
            net_Na_msd_bulk = np.array((0.0, 0.0, 0.0))
            Na_msd_com_bulk = np.array((0.0, 0.0, 0.0))
            net_Na_msd_com_bulk = np.array((0.0, 0.0, 0.0))

            for frame in range(i, n_frames): 
                #print(frame)
                data = pipeline.compute(frame)
                Na_msd_bulk += Calculate_MSD_grain_raw(data, Na_type_id)
                net_Na_msd_bulk += Calculate_netMSD_grain_raw(data, Na_type_id)
                Na_msd_com_bulk += Calculate_MSD_grain_com(data, Na_type_id)
                net_Na_msd_com_bulk += Calculate_netMSD_grain_com(data, Na_type_id)
                
                #print(Na_msd_bulk)

            Na_MSD_bulk.append(Na_msd_bulk/(n_frames - i))
            net_Na_MSD_bulk.append(net_Na_msd_bulk/(n_frames -i ))
            Na_MSD_com_bulk.append(Na_msd_com_bulk/(n_frames - i))
            net_Na_MSD_com_bulk.append(net_Na_msd_com_bulk/(n_frames -i ))

        #print(Na_MSD_bulk[0])


        #Export results to txt file
        t = np.arange(len(Na_MSD_bulk), dtype = float)
        t*=(dt*every_nth_frame)
        np.savetxt(file_2, np.column_stack( (t, Na_MSD_bulk[:][0], Na_MSD_bulk[:][1], Na_MSD_bulk[:][2], 
                                   net_Na_MSD_bulk[:][0], net_Na_MSD_bulk[:][1], net_Na_MSD_bulk[:][2], 
                                   Na_MSD_com_bulk[:][0], Na_MSD_com_bulk[:][1], Na_MSD_com_bulk[:][2], 
                                   net_Na_MSD_com_bulk[:][0], net_Na_MSD_com_bulk[:][1], net_Na_MSD_com_bulk[:][2],
                                  )), delimiter = " ", header = "Timestep(ps),MSDx_bulk_raw(A^2),MSDy_bulk_raw(A^2),MSDz_bulk_raw(A^2),netMSDx_bulk_raw(A^2),netMSDy_bulk_raw(A^2),netMSDz_bulk_raw(A^2),MSDx_bulk_com(A^2),MSDy_bulk_com(A^2),MSDz_bulk_com(A^2),netMSDx_bulk_com(A^2),netMSDy_bulk_com(A^2),netMSDz_bulk_com(A^2)")
