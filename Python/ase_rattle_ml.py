from statistics import stdev
import numpy as np

# Function that distort the cell, rattle atoms etc.  
# Anisotropic deformation 

def deform_aniso(atoms, seed=None, pos_shake=0.05, cell_shake=0.05):
    
    if seed is not None:
        np.random.seed(seed)
    #seed = np.random.randint()
    seed = np.random.seed(100)
    atoms.rattle(stdev=0.05, seed = seed)                      # Randomly displace the atoms 
    orig_pos = atoms.get_positions()                           #  
    orig_cell = atoms.get_cell()                               # 
    dpos = np.random.randn(*orig_pos.shape) * pos_shake        # Shaking the position of atoms 
    dcell_matrix = np.eye(3)+np.random.randn(3,3)*cell_shake   # Shake matrix to dicplace volume 
    new_cell = np.dot(orig_cell, dcell_matrix)
    new_pos = orig_pos+dpos
    new_atoms = atoms.copy()
    new_atoms.set_positions(new_pos)
    new_atoms.set_cell(new_cell,scale_atoms=True)
    return new_atoms

def deform_iso(atoms, seed=None, pos_shake=0.05, cell_shake=0.05):
    
    if seed is not None:
        np.random.seed(seed)
    #seed = np.random.randint()
    seed = np.random.seed(100)
    atoms.rattle(stdev=0.05, seed = seed)                       # Randomly displace the atoms 
    orig_pos = atoms.get_positions()                            #  
    orig_cell = atoms.get_cell()                                # 
    dpos = np.random.randn(*orig_pos.shape) * pos_shake         # Shaking the position of atoms 
    #dcell_matrix = np.eye(3)+np.random.randn(3,3)*cell_shake   # Shake matrix to displace volume anisotropicaly 
    dcell_matrix = np.eye(3)*cell_shake                         # Shake matrix to displace volume isotropicaly   
    new_cell = np.dot(orig_cell, dcell_matrix)
    new_pos = orig_pos+dpos
    new_atoms = atoms.copy()
    new_atoms.set_positions(new_pos)
    new_atoms.set_cell(new_cell,scale_atoms=True)
    return new_atoms


def rattle(atoms, seed=None):
    
    if seed is not None:
        np.random.seed(seed)
    #seed = np.random.randint()
    seed = np.random.seed(100)
    new_atoms = atoms.rattle(stdev=0.05, seed = seed)                       # Randomly displace the atoms 
    return new_atoms
