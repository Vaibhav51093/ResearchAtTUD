from ovito.io import import_file
from ovito.modifiers import SelectTypeModifier
from ovito.modifiers import SpatialBinningModifier
from ovito.io import export_file
import numpy as np
import operator

### Parameters ########################################################################################################################################
in_filename       = 'XDATCAR'
out_filename      = 'Li_density.xyz' # There will be scaled RDFs (all partial sum up to the total RDF) and unscaled (as if only those particles were in the box)
EquilibrationTime = 500   # Skip this number of STEPS from the XDATCAR. Not necessary fs -> remember NBLOCK from INCAR!
bins_x            = 192    # number of bins in x direction
bins_y            = 192    # number of bins in y direction
bins_z            = 256   # number of bins in z direction
AtomType          = 'Li'  # Species to be binned
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


### Set up the binning operation. Binning can be done over "Selection" to count if any Li (is selected, hence: Selection = 1) is in the bin or not.
binning     = SpatialBinningModifier(
        bin_count_x   = bins_x,
        bin_count_y   = bins_y,
        bin_count_z   = bins_z,
        property      = "Selection",
        only_selected = True,
        direction     = SpatialBinningModifier.Direction.XYZ,
        reduction_operation = SpatialBinningModifier.Operation.SumVol
        )
pipeline.modifiers.append(binning)


### Compute pipeline
print("Skipping first {} frames of {}".format(EquilibrationTime,pipeline.source.num_frames))
grid = np.zeros(bins_x * bins_y * bins_z ,)

for frame in range(EquilibrationTime,pipeline.source.num_frames):
#for frame in range(EquilibrationTime,EquilibrationTime+1000):
    if frame % 1000 == 0:
        print("Computing frame {} of {}".format(frame,pipeline.source.num_frames),flush=True)
    # Evaluate pipeline to let the modifier compute the RDF of the current frame:
    data = pipeline.compute(frame)
    #print(data.grids)
    Li_selection = data.grids['binning[Selection]']['Selection']
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
    
    lattice_string=""
    for i in a:
        lattice_string=lattice_string + "{:.8f} ".format(i)
    for i in b:
        lattice_string=lattice_string + "{:.8f} ".format(i)
    for i in c:
        lattice_string=lattice_string + "{:.8f} ".format(i)
    f.write('Lattice="{}" Properties=pos:R:3:Density:R:1 Time=0.0\n'.format(lattice_string.rstrip()))
    
    for z in range(bins_z):
        for y in range(bins_y):
            for x in range(bins_x):
                Coordinates = np.array([0,0,0])
                Coordinates =  Coordinates + x*a/bins_x +  y*b/bins_y + z*c/bins_z
                f.write( "{:.8f}\t{:.8f}\t{:.8f}\t{:.8f}\n".format(Coordinates[0],Coordinates[1],Coordinates[2],grid[z][y][x]) )

