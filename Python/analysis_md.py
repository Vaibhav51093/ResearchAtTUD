from ovito.modifiers import CalculateDisplacementsModifier
from ovito.modifiers import CoordinationNumberModifier
from ovito.modifiers import VoronoiAnalysisModifier
from ovito.modifiers import PythonScriptModifier
from ovito.io import import_file
#from scipy import stats as st
import os
import numpy as np 


class ovito_tool:

    def cal_msd(step=1.0, dump_write=100, step_ovito=10, file = 'cryt_msd_773k.txt', path='/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass/boilot_623K_2_05_cryt_773_k_ps_new_hdf5/boilot_623K_2_05_cryt_773_k_ps_new/dump_nvt_prod.out'):
        # Simulation Parameters
        POTIM = step                    # Time step in fs
        NBLOCK = dump_write             # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?

        # MSD Evaluation
        Stepsize = step_ovito

        pipeline = import_file(path)
        # Print the list of input particle type

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

        # Extract Atom Types
        AtomTypes = []
        data = pipeline.compute()
        for Atom in data.particles['Particle Type']:
            if data.particles['Particle Type'].type_by_id(Atom).name not in AtomTypes:
                AtomTypes.append(data.particles['Particle Type'].type_by_id(Atom).name)


        # Extract number of Atoms for each Atom Type
        AtomTypesCount = []
        for Type in AtomTypes:
            count = 0
            for Atom in data.particles['Particle Type']:
                if data.particles['Particle Type'].type_by_id(Atom).name == Type:
                    count += 1
            AtomTypesCount.append(count)


        # Set up Displacement Modifier and append to pipeline
        displ_mod    = CalculateDisplacementsModifier(
            use_frame_offset  = True,
            frame_offset      = -Stepsize)
        pipeline.modifiers.append(displ_mod)


        # Set up empty array to store intermediate displacements
        Displacements = []
        for Atom in data.particles['Particle Type']:
            Displacements.append([0.0,0.0,0.0])


        # Extract Displacements
        Time_and_MSD = [] #Structure will be:  [ Time, Li_Total, [Li_X,Li_Y,Li_Z], P_Total, [P_X,P_Y,P_Z], ... ]
        print("Starting MSD calculation for complete Trajectory")
        for frame in range(Stepsize,pipeline.source.num_frames,Stepsize):

            # Just to see a bit of update
            if frame % 1000 == 0:
                print("    Computing frame %d of %d" % (frame,pipeline.source.num_frames),flush=True)


            # Set up temporary List that is added to Time_and_MSD list after each evaluation
            # Structure is:  [ Time, Li_Total, [Li_X,Li_Y,Li_Z], P_Total, [P_X,P_Y,P_Z], ... ]
            Time_and_MSD_temporary = [frame*NBLOCK*POTIM/1000]
            for Atom in AtomTypes:
                Time_and_MSD_temporary.append(0.0)
                Time_and_MSD_temporary.append([0.0,0.0,0.0])


            # Evaluate pipeline to let the modifier compute the Displacements of the current frame:
            data = pipeline.compute(frame)

            # Each iteration: Update the displacement vector of each atom, check its atom type and add the MSD to the Time_and_MSD_temporary list.
            for (AtomID, Atom, displ) in zip(range(0, sum(AtomTypesCount)), data.particles['Particle Type'], data.particles['Displacement']):
                Displacements[AtomID] = np.add(Displacements[AtomID],displ)
                for ListPosition,Type in enumerate(AtomTypes):
                    if data.particles['Particle Type'].type_by_id(Atom).name == Type:
                        Time_and_MSD_temporary[1+ListPosition*2]      += (Displacements[AtomID][0]**2 + Displacements[AtomID][1]**2 + Displacements[AtomID][2]**2) / AtomTypesCount[ListPosition]
                        Time_and_MSD_temporary[1+ListPosition*2+1][0] += (Displacements[AtomID][0]**2                                                            ) / AtomTypesCount[ListPosition]
                        Time_and_MSD_temporary[1+ListPosition*2+1][1] += (                              Displacements[AtomID][1]**2                              ) / AtomTypesCount[ListPosition]
                        Time_and_MSD_temporary[1+ListPosition*2+1][2] += (                                                            Displacements[AtomID][2]**2) / AtomTypesCount[ListPosition]
            Time_and_MSD.append(Time_and_MSD_temporary)


        # Header File
        header = '#MSD evaluated over all available data. Time in ps and MSDs in Ang^2\n#Time'
        for Type in AtomTypes:
            header = header + '\t  ' + Type + '\t  ' + Type + '_X\t  ' + Type + '_Y\t  ' + Type + '_Z'
        header = header + '\n  0.000\t'
        for Type in AtomTypes:
            header = header + '  0.000\t' + '  0.000\t' + '  0.000\t' + '  0.000\t'
        header = header + '\n'

        with open(file, 'w') as f:
            f.write(header)
            for DataPoints in Time_and_MSD:
                f.write('{:7.3f}'.format(DataPoints[0])) # Time
                for ListPosition in range(0,len(AtomTypes)):
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2  ]   ))
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2+1][0]))
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2+1][1]))
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2+1][2]))
                f.write('\n')

