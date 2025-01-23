from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

# Provide the lattice parameters (a, b, c) in Angstroms.
lattice_parameters = [2.8815, 2.8815, 11.244]

# Fixed lattice angles (alpha, beta, gamma) in degrees.
alpha, beta, gamma = 90.0, 90.0, 120.0

# Provide the list of atomic symbols in the same order as the atomic positions.
atomic_symbols = ['Na', 'Na', 'Mn', 'Mg', 'O']

# Provide the atomic positions in fractional coordinates.
# Make sure the fractional coordinates are between 0 and 1.
atomic_positions = [(0.000, 0.000, 0.25), (0.333, 0.6667, 0.7500), (0.000, 0.000, 0.000), (0.000, 0.000, 0.000), (0.333, 0.6667, 0.0901)]

# Provide the occupancies for each site (corresponding to each atomic position).
# Occupancies for Mn and Mg should add up to less than 1.0 for each site.
occupancies = [0.0152, 0.0277, 0.0770*8, 0.0082*7, 0.1666]

# Create the lattice matrix manually with fixed angles.
lattice_matrix = Lattice.from_parameters(a=lattice_parameters[0], b=lattice_parameters[1], c=lattice_parameters[2],
                                         alpha=alpha, beta=beta, gamma=gamma)

# Create the Structure object
structure = Structure(lattice_matrix, atomic_symbols, atomic_positions, occupancies)

# Get the primitive cell structure
primitive_structure = structure.get_primitive_structure()

# Write the CIF file for the primitive cell
primitive_structure.to(filename="primitive_crystal_structure.cif")

