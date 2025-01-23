import time
start_time = time.time()

from random import random
from ase.io import read, write
from ase.db import connect
from ase.optimize import GPMin
from gpaw import GPAW, FermiDirac
from ase.optimize import BFGS
from ase.ga.data import DataConnection
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.offspring_creator import OperationSelector
from ase.ga.standardmutations import MirrorMutation
from ase.ga.standardmutations import RattleMutation
from ase.ga.standardmutations import PermutationMutation
import argparse

parser = argparse.ArgumentParser()
defval = 0.3
parser.add_argument('-m', '--mutation_probability', type=float, default=defval,
                    help='Mutation probability (default: {:.2f})'.format(defval))
defval = 80
parser.add_argument('-n', '--ncandidates_to_test', type=int, default=defval,
                    help='Number of candidates to test (default: {})'.format(defval))
args = parser.parse_args()

# Initialize the different components of the GA
#da = DataConnection('gadb.db')
#da = read('christmas-tree_last.xyz', format='xyz') 
# Determine size of starting population
# Read ase.io 
#db = connect('gadb.db').get('id=21')
db = read('half-decahedron_second_last.xyz', format='extxyz')
#population_size = db.count('natoms>0')
#population_size = len(db.get_positions())
#print('Running with population_size = {}'.format(population_size))




# Define the calculator
calc = GPAW(nbands=10,
            h=0.25,
            txt='out.txt',
            occupations=FermiDirac(0.05),
            setups={'Na': '1'},
            mode='fd',
            xc='PBE')
            #basis='dzp')

# Relax the new candidate
db.set_calculator(calc)
dyn = GPMin(db, trajectory='relax_new.traj', logfile='relax_new.log')
dyn.run(fmax=0.02, steps=100)
#db.info['key_value_pairs']['raw_score'] = -db.get_potential_energy()
#db.add_relaxed_step(db)

#write('traj.traj', db.get_all_relaxed_candidates())
#opt = BFGS(da)
#write('all_candidates.traj', da)
#da.get_potential_energy()
print("--- {:1f} seconds ---".format(time.time() - start_time))
