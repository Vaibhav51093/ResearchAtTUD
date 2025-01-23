from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
import operator
import numpy as np
from scipy import integrate
from scipy.constants import codata

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev

def na_density_na(path_line = "/nfshome/deshmukh/vaibhav/schulze_project/gb/IDp02/00-minimize/01-heat-773/02-diffusion/structure.dump.*", id=1, bin_x=700, x=False, y=False):

    pipeline = import_file(path_line)

    def setup_particle_types(frame, data):
                types = data.particles_.particle_types_
                types.type_by_id_(1).name = "Zr"
                types.type_by_id_(2).name = "Si"
                types.type_by_id_(3).name = "P"
                types.type_by_id_(4).name = "Na"
                types.type_by_id_(5).name = "O"
                types.type_by_id_(5).radius = 0.3
        
    pipeline.modifiers.append(setup_particle_types)

    # Wrap at periodic boundaries:
    pipeline.modifiers.append(WrapPeriodicImagesModifier())

    #if type == True:
    id = id

    pipeline.modifiers.append(SelectTypeModifier(types = {id}))

    # Spatial binning:

    if x == True:
        pipeline.modifiers.append(SpatialBinningModifier(
            property = 'Selection', 
            reduction_operation = SpatialBinningModifier.Operation.SumVol, 
            direction = SpatialBinningModifier.Direction.X, 
            bin_count = (bin_x, 80, 46)))
        data = pipeline.compute()
        cell = data.cell[:,0][0]
    elif y == True:
        pipeline.modifiers.append(SpatialBinningModifier(
            property = 'Selection', 
            reduction_operation = SpatialBinningModifier.Operation.SumVol, 
            direction = SpatialBinningModifier.Direction.Y, 
            bin_count = (bin_x, 80, 46)))
        data = pipeline.compute()
        cell = data.cell[:,1][1]
    else:
        pipeline.modifiers.append(SpatialBinningModifier(
            property = 'Selection', 
            reduction_operation = SpatialBinningModifier.Operation.SumVol, 
            direction = SpatialBinningModifier.Direction.Z, 
            bin_count = (bin_x, 80, 46)))
        data = pipeline.compute()
        cell = data.cell[:,2][2]
        
    
    # Time averaging:
    #pipeline.modifiers.append(TimeAveragingModifier(operate_on = 'table:binning'))

    charge_density = []

    for i in range(pipeline.source.num_frames):

        # Calculate electrostatic potential:
        data = pipeline.compute(frame=i)
        charge_density.append(data.tables['binning']['Selection'])

    charge_density = np.average(charge_density, axis=0)
    
    return charge_density, np.linspace(0,cell,bin_x)

def compute_efield_esp(charge_density, data):
    
    '''This function computes the electric field and electrostatic potential from the charge density
       using the trapezoidal rule for integration. The function returns the electric field and electrostatic potential
       as well as the spatial grid on which the charge density is defined. The function also shifts the electric field
       and potential to have zero mean. The function takes the charge density and the spatial grid as input and returns'''
       
    
    e_filed_tt = 14.3997584 * integrate.cumtrapz(charge_density, data, initial=0) # The scale factor of 14.3997584 is used to convert the electric field from atomic units (a.u.) to volts per angstrom (V/Å).
    e_filed_tt = e_filed_tt - np.mean(e_filed_tt)  # Shift the electric field to zero mean

    potential_tt = -integrate.cumtrapz(e_filed_tt, data, initial=0)
    potential_tt = potential_tt - np.mean(potential_tt) # Shift the potential to zero mean

    return e_filed_tt, potential_tt, data