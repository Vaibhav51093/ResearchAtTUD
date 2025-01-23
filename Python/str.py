from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

## Provide the lattice parameters (a, b, c, alpha, beta, gamma) in Angstroms and degrees.
#lattice_parameters = [[2.8815, 0, 0], [0, 2.882, 0], [0, 0, 11.244]]

# Provide the lattice parameters (a, b, c) in Angstroms.
lattice_parameters = [2.8815, 2.882, 11.244]

# Fixed lattice angles (alpha, beta, gamma) in degrees.
alpha, beta, gamma = 90.0, 90.0, 120.0



# Provide the list of atomic symbols in the same order as the atomic positions.
atomic_symbols = ['Na', 'Na', 'Mn', 'Mg', 'O']

# Provide the atomic positions in fractional coordinates.
# Make sure the fractional coordinates are between 0 and 1.
atomic_positions = [(0.000, 0.000, 0.25), (0.333, 0.6667, 0.7500), (0.000, 0.000, 0.000), (0.000, 0.000, 0.000), (0.333, 0.6667, 0.0901)]

# Provide the occupancies for each site (corresponding to each atomic position).
occupancies = [0.0152, 0.0277, 0.0770, 0.0082, 0.1666]

lattice_matrix = Lattice.from_parameters(a=lattice_parameters[0], b=lattice_parameters[1], c=lattice_parameters[2],
                                         alpha=alpha, beta=beta, gamma=gamma)


# Create the Structure object
#structure = Structure(Lattice(lattice_parameters), atomic_symbols, atomic_positions, occupancies)


# Create the Structure object
structure = Structure(lattice_matrix, atomic_symbols, atomic_positions, occupancies)


# Write the CIF file
structure.to(filename="crystal_structure.cif")

