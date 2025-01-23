from random import random
from ase.io import read, write
from ase.ga.data import PrepareDB

atom_numbers = 6

da = read('half-decahedron_second_last.xyz', format='xyz')
PrepareDB(db_file_name='half-decahedron_second_last.db',
                  simulation_cell=da,
                  stoichiometry=atom_numbers)