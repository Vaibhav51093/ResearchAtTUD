# Bond analysis script OVITO
from ovito.data import *
import numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import *
import itertools

# Import lammps dump file  
pipeline = import_file("dump_amo_723.out", multiple_frames = True)

# Manual modifications of the imported data objects:
def modify_pipeline_input(frame: int, data: DataCollection):
    data.particles_.particle_types_.type_by_id_(1).color = (0.6666666865348816, 0.0, 1.0)
    data.particles_.particle_types_.type_by_id_(1).radius = 0.2
    data.particles_.particle_types_.type_by_id_(2).color = (1.0, 1.0, 0.0)
    data.particles_.particle_types_.type_by_id_(2).radius = 0.5
    data.particles_.particle_types_.type_by_id_(3).radius = 0.2
    data.particles_.particle_types_.type_by_id_(4).radius = 0.2
    data.particles_.particle_types_.type_by_id_(5).radius = 0.2
pipeline.modifiers.append(modify_pipeline_input)


# Asign the names 
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

# Create bonds:
mod = CreateBondsModifier()
mod.mode = CreateBondsModifier.Mode.Pairwise
pipeline.modifiers.append(mod)
mod.set_pairwise_cutoff('P', 'O', 2.0)        # Official bond length in Padone for p-o 
mod.set_pairwise_cutoff('Si', 'O', 2.2)       # -----------------""-------------------
mod.set_pairwise_cutoff('Zr', 'O', 2.5)       # -----------------""-------------------

for frame in range(pipeline.source.num_frames):
    data = pipeline.compute(frame)
    bonds_topology = np.array(data.particles.bonds.topology)
    # Bond enumerator 
    bonds_enum = BondsEnumerator(data.particles.bonds)
    #Request the particle type property
    particle_types = data.particles["Particle Type"]
    ptypes = data.particles.particle_type

    possible_combinations_of_bond_pairs= [np.array((bond1, bond2)) for bond1,bond2 in itertools.combinations(bonds_topology, r=2) if (bond1[0] in bond2) or (bond1[1] in bond2)]
    possible_combinations_of_particle_type_triplets = [ptypes[np.unique(topo)] for topo in possible_combinations_of_bond_pairs]

    triplet_type_ids, count = np.unique(possible_combinations_of_particle_type_triplets, axis = 0, return_counts = True)
    #print(triplet_type_ids, count)
    triplet_type_names = [(ptypes.type_by_id(t[0]).name, ptypes.type_by_id(t[1]).name, ptypes.type_by_id(t[2]).name) for t in triplet_type_ids  ]
    #print(triplet_type_names, count)
    count_ZrO2 = []
    count_SiP2 = []
    for i in triplet_type_names:
        if frame == 5001:
            if i == ('O', 'O', 'Zr') or ('O', 'Zr', 'O') or ('Zr', 'O', 'O'):
                count_ZrO2.append(i)
                print(i)
        else:
            continue
    avg_count = len(count_ZrO2)/len(range(pipeline.source.num_frames))