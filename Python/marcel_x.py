from fileinput import filename
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
#from scipy import stats as st
import tarfile
import lammps_logfile
import os
import codecs
utf8reader = codecs.getreader('utf-8')
import os
import numpy as np 
from scipy import integrate
from scipy.constants import codata
import pandas as pd
from ase.io import read, write
from scipy.optimize import curve_fit
from pyiron import Project
from pyiron import ase_to_pyiron, pyiron_to_ase
import shutil
import glob
from scipy import stats
from scipy.optimize import curve_fit
import numpy as np 
import scipy.constants as const
import sys

sys.path.append(os.path.abspath("/nfshome/deshmukh/vaibhav/scripts"))
import analysis_msd as m

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev

def diff_3(path_line='/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_ordered_glass/boilot_623K_2_05_cryt_523_k_ps_new_hdf5/boilot_623K_2_05_cryt_523_k_ps_new/dump_nvt_prod.out',file_2='msd_na_y_45_jochen.txt',step=1.0, dump_write=100, step_ovito=10):

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

    # Rotation of the box to 30 degree clockwise to extraxt Dyy 
    # Create a rotation matrix
    #theta = 30.0 # Angle of rotation in degrees
    #theta = np.radians(theta) # Convert to radians
    #R = np.array([[1, 0, 0, 0],
    #              [0, np.cos(theta), -np.sin(theta), 0],
    #              [0, np.sin(theta), np.cos(theta), 0]])#

    # Apply the rotation matrix to the simulation box

    modifier = AffineTransformationModifier(transformation = [[0.7071067811865476, 0.0, 0.7071067811865475, -6.176977571936642], [0.0, 1.0, 0.0, 0.0], [-0.7071067811865475, 0.0, 0.7071067811865476, 8.205565304527283]])
    #pipeline.modifiers.append(AffineTransformationModifier(transformation = [[0.8660254037844387, 0.49999999999999994, 0.0, 0.0], [-0.49999999999999994, 0.8660254037844387, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]))
    # Displacement vectors

    pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

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

    def modify3(frame: int, data: DataCollection, type_id_Na = 2):

        ptype = data.particles.particle_type

        time = frame*NBLOCK*POTIM/1000

        bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
        bulk_displ = data.particles['Displacement'][(ptype == type_id_Na)]
        
        msd_x, msd_y, msd_z = np.mean(bulk_displ_rel**2, axis = 0)
        data.attributes["MSDx_bulk_com"] = msd_x  # com = center of mass /displacements are corrected for com shift of host matrix
        data.attributes["MSDy_bulk_com"] = msd_y  # MSD = tracer diffucivity
        data.attributes["MSDz_bulk_com"] = msd_z
        data.attributes["MSDav_bulk_com"] = msd_x+msd_y+msd_z

        net_x, net_y, net_z = np.sum(bulk_displ_rel, axis = 0)
        data.attributes["netMSDx_bulk_com"] = net_x**2 / len(bulk_displ_rel) # netMSD = ionic diffusivity
        data.attributes["netMSDy_bulk_com"] = net_y**2 / len(bulk_displ_rel)
        data.attributes["netMSDz_bulk_com"] = net_z**2 / len(bulk_displ_rel)
        data.attributes["netMSDav_bulk_com"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ_rel)

        msd_x, msd_y, msd_z = np.mean(bulk_displ**2, axis = 0)
        data.attributes["MSDx_bulk_raw"] = msd_x # raw = raw data/shift of host matrix is not considered
        data.attributes["MSDy_bulk_raw"] = msd_y 
        data.attributes["MSDz_bulk_raw"] = msd_z
        data.attributes["MSDav_bulk_raw"] = msd_x+msd_y+msd_z

        net_x, net_y, net_z = np.sum(bulk_displ, axis = 0)
        data.attributes["netMSDx_bulk_raw"] = net_x**2 / len(bulk_displ)
        data.attributes["netMSDy_bulk_raw"] = net_y**2 / len(bulk_displ)
        data.attributes["netMSDz_bulk_raw"] = net_z**2 / len(bulk_displ)
        data.attributes["netMSDav_bulk_raw"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ)

        data.attributes["Timestep"] = time
   
        
    pipeline.modifiers.append(modify3)

    export_file(pipeline, file_2, "txt/attr",
    columns=["Timestep", "MSDx_bulk_raw", "MSDy_bulk_raw", "MSDz_bulk_raw","MSDav_bulk_raw","netMSDx_bulk_raw", "netMSDy_bulk_raw", "netMSDz_bulk_raw","netMSDav_bulk_raw","MSDx_bulk_com", "MSDy_bulk_com", "MSDz_bulk_com", "MSDav_bulk_com","netMSDx_bulk_com", "netMSDy_bulk_com", "netMSDz_bulk_com","netMSDav_bulk_com"],
    multiple_frames=True)