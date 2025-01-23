from ovito.data import *
import numpy as np
def modify(frame, data):
   
    #Get the particle type property
    particle_types = data.particles.particle_types
    #Get bonds topology and particle positions to compute all bond vectors present in your system
    topology = data.particles.bonds.topology
    positions = data.particles.positions
    bond_vectors = positions[topology[:,1]] - positions[topology[:,0]]
    bond_vectors += np.dot(data.cell[:3,:3], data.particles.bonds.pbc_vectors.T).T
    #Store this information as custom bond property, so they appear in the Data inspector
    data.particles_.bonds_.create_property("bond_vectors", data = bond_vectors)    
    # Set up the bonds enumerator object.
    bonds_enum = BondsEnumerator(data.particles.bonds)
    # Counter
    bond_count = 0
 
    # Loop over atoms.   
    for particle_index in range(data.particles.count):
        # If the particle is a Hw atom continue and loop over bonds of current atom.
        if particle_types[particle_index] == 5:
            #Make sure that the current particle is connected to 2 bonds 
            bond_indices_list = [bond_index for bond_index in bonds_enum.bonds_of_particle(particle_index)]
            if len(bond_indices_list) == 2:
                                    
                neighbor_particle_indices = []
                for bond_index in bonds_enum.bonds_of_particle(particle_index):
                    # Obtain the indices of the two particles connected by the bond:
                    a = topology[bond_index, 0]
                    b = topology[bond_index, 1]
                    if a == particle_index:
                        neighbor_particle_indices.append(b)
                    else:
                        neighbor_particle_indices.append(a) 
                    
                #Verify that the particle types of the 2 neighboring atoms are actually to type Ow
           
                if( (particle_types[neighbor_particle_indices] == np.array([4,4])).all() ):
                    #Calculate angle
                    v_1 = bond_vectors[bond_indices_list[0]] 
                    v_2 = bond_vectors[bond_indices_list[1]]
                    angle = np.arccos(np.clip(np.dot(v_1, v_2)/np.linalg.norm(v_1)/np.linalg.norm(v_2), -1.0,1.0 ))
                   
                    if np.degrees(angle) >= 120  :
                        bond_count+=1
    print("# Ow-Hw-Ow bonds:{}".format(bond_count))
    #Store result as global attribute so it appears in the Data inspector in the GUI
    data.attributes["Ow-Hw-Ow"] = bond_count