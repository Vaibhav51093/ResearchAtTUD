from ovito.io import import_file
from ovito.modifiers import SpatialBinningModifier
import numpy as np

### Parameters
in_filename = 'structure.dump.*'
out_filename = 'charge_density.xyz'
EquilibrationTime = 0
bins_x = 50
bins_y = 50
bins_z = 300

### Load a simulation trajectory consisting of several frames:
pipeline = import_file(in_filename)

# Spatial binning:
pipeline.modifiers.append(SpatialBinningModifier(
    property='Charge',
    direction=SpatialBinningModifier.Direction.XYZ,
    reduction_operation=SpatialBinningModifier.Operation.SumVol,
    bin_count=(bins_x, bins_y, bins_z)))

### Compute pipeline
print("Skipping first {} frames of {}".format(EquilibrationTime, pipeline.source.num_frames))
grid = np.zeros((bins_z, bins_y, bins_x))

for frame in range(EquilibrationTime, pipeline.source.num_frames):
    if frame % 100 == 0:
        print("Computing frame {} of {}".format(frame, pipeline.source.num_frames), flush=True)
    data = pipeline.compute(frame)
    charge_density = data.grids['binning'] # Access the array data of the VoxelGrid object
    #charge_density = charge_density#.astype('float64')  # Ensure the charge density array is of the correct type
    charge_density = charge_density['Charge'][:].reshape((bins_z, bins_y, bins_x))
    grid += charge_density

    #grid.append(charge_density)

grid /= (pipeline.source.num_frames - EquilibrationTime)

a = np.array([data.cell[0][0], data.cell[1][0], data.cell[2][0]])
b = np.array([data.cell[0][1], data.cell[1][1], data.cell[2][1]])
c = np.array([data.cell[0][2], data.cell[1][2], data.cell[2][2]])

with open(out_filename, 'w') as f:
    f.write("{}\n".format(bins_x * bins_y * bins_z))
    f.write("X\tY\tZ\tCharge_density\n")
    for z in range(bins_z):
        for y in range(bins_y):
            for x in range(bins_x):
                Coordinates = np.array([0, 0, 0])
                Coordinates = Coordinates + x * a / bins_x + y * b / bins_y + z * c / bins_z
                f.write("{:.4f}\t{:.4f}\t{:.4f}\t{:.8f}\n".format(Coordinates[0], Coordinates[1], Coordinates[2], grid[z][y][x]))
