from ovito.io import import_file
from ovito.modifiers import SelectTypeModifier
from ovito.modifiers import SpatialBinningModifier
from ovito.modifiers import ComputePropertyModifier
from ovito.modifiers import WrapPeriodicImagesModifier
from ovito.io import export_file
import numpy as np
import operator

### Parameters ########################################################################################################################################
in_filename       = '/nfshome/deshmukh/vaibhav/shanka_cluster/deshmukh2/nas_0/623/dump_nvt_prod_3.out'
out_filename      = 'Na_density.xyz' # There will be scaled RDFs (all partial sum up to the total RDF) and unscaled (as if only those particles were in the box)
EquilibrationTime = 0     # Skip this number of STEPS from the XDATCAR. Not necessary fs -> remember NBLOCK from INCAR!
bins_x            = 200   # number of bins in x direction
bins_y            = 200   # number of bins in y direction
bins_z            = 200   # number of bins in z direction
AtomType          = 2     # Species to be binned, sodium 
########################################################################################################################################################


### Load a simulation trajectory consisting of several frames:
pipeline = import_file(in_filename)

# Manual modifications of the imported data objects:
def setup_particle_types(frame, data):
    types = data.particles_.particle_types_
    types.type_by_id_(1).name = "O"
    types.type_by_id_(2).name = "Na"
    types.type_by_id_(3).name = "P"
    #types.type_by_id_(4).name = "Si"
    types.type_by_id_(5).name = "Zr"
        
pipeline.modifiers.append(setup_particle_types)

def modify1(frame, data):

    charge_dict = {"Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
    mass_dict = {"Zr": 91.224,
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

## Wrap at periodic boundaries:
#pipeline.modifiers.append(WrapPeriodicImagesModifier())

### Select Atom Type and add to pipeline
sel_modifier = SelectTypeModifier(
        operate_on = "particles",
        property   = "Particle Type",
        )
sel_modifier.types = {AtomType}
pipeline.modifiers.append(sel_modifier)
pipeline.modifiers.append(WrapPeriodicImagesModifier())

# Compute property:
pipeline.modifiers.append(ComputePropertyModifier(
    expressions = ('1',), 
    output_property = 'Charge', 
    only_selected = False))

### Set up the binning operation. Binning can be done over "Selection" to count if any Li (is selected, hence: Selection = 1) is in the bin or not.
binning     = SpatialBinningModifier(
        bin_count = (bins_x, bins_y, bins_z),
        #bin_count_x   = bins_x,
        #bin_count_y   = bins_y,
        #bin_count_z   = bins_z,
        property      = "Selection",
        only_selected = False,
        direction     = SpatialBinningModifier.Direction.XYZ,
        reduction_operation = SpatialBinningModifier.Operation.Sum
        )
pipeline.modifiers.append(binning)

### Compute pipeline
print("Skipping first {} frames of {}".format(EquilibrationTime,pipeline.source.num_frames))
grid = np.zeros(bins_x * bins_y * bins_z ,)

for frame in range(EquilibrationTime,pipeline.source.num_frames):
#for frame in range(EquilibrationTime,EquilibrationTime+1000):
    if frame % 100 == 0:
        print("Computing frame {} of {}".format(frame,pipeline.source.num_frames),flush=True)
    # Evaluate pipeline to let the modifier compute the RDF of the current frame:
    data = pipeline.compute(frame)
    #print(data.grids)
    Li_selection = data.grids['binning']['Selection']
    grid+=Li_selection
    
grid/= (pipeline.source.num_frames-EquilibrationTime)
grid = np.reshape( grid, (bins_z,bins_y,bins_x) )

a = np.array([data.cell[0][0],data.cell[1][0],data.cell[2][0]])
b = np.array([data.cell[0][1],data.cell[1][1],data.cell[2][1]])
c = np.array([data.cell[0][2],data.cell[1][2],data.cell[2][2]])

print(a)
print(b)
print(c)

with open(out_filename, 'w') as f:
    f.write("{}\n".format(bins_x*bins_y*bins_z))
    f.write("X\tY\tZ\tLi_density\n")
    for z in range(bins_z):
        for y in range(bins_y):
            for x in range(bins_x):
                Coordinates = np.array([0,0,0])
                Coordinates =  Coordinates + x*a/bins_x +  y*b/bins_y + z*c/bins_z
                f.write( "{:.4f}\t{:.4f}\t{:.4f}\t{:.8f}\n".format(Coordinates[0],Coordinates[1],Coordinates[2],grid[z][y][x]) )
