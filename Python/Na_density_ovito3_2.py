from ovito.io import import_file
from ovito.modifiers import SelectTypeModifier
from ovito.modifiers import SpatialBinningModifier
from ovito.modifiers import ComputePropertyModifier
from ovito.modifiers import WrapPeriodicImagesModifier
from ovito.io import export_file
import numpy as np
import operator

### Parameters ########################################################################################################################################
in_filename       = '/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_2_2/minimization/hena_2_2/hena_1_struct_eq_big_623k_6/dump_nvt_prod.out'
out_filename      = 'na_big_cryt_623k_6.xyz' # There will be scaled RDFs (all partial sum up to the total RDF) and unscaled (as if only those particles were in the box)
EquilibrationTime = 0     # Skip this number of STEPS from the XDATCAR. Not necessary fs -> remember NBLOCK from INCAR!
bins_x            = 50    # number of bins in x direction
bins_y            = 50    # number of bins in y direction
bins_z            = 50   # number of bins in z direction
AtomType          = 2     # Species to be binned
########################################################################################################################################################


### Load a simulation trajectory consisting of several frames:
pipeline = import_file(in_filename)


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
    output_property = 'Selection', 
    only_selected = True))

### Set up the binning operation. Binning can be done over "Selection" to count if any Li (is selected, hence: Selection = 1) is in the bin or not.
binning     = SpatialBinningModifier(
        bin_count = (bins_x, bins_y, bins_z),
        #bin_count_x   = bins_x,
        #bin_count_y   = bins_y,
        #bin_count_z   = bins_z,
        property      = "Selection",
        only_selected = True,
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

