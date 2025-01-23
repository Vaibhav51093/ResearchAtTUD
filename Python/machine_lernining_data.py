from pyiron import Project, ase_to_pyiron
import matplotlib.pyplot as plt
import numpy as np
from pyiron import Project
from ase.io import read, write
from pyiron import ase_to_pyiron
import ase
import os
import time



class shake_rattle_():

    # Function to change the atoms positions and others 
    def shake_atoms(atoms, seed=None, pos_shake=0.05, cell_shake=0.05):
        
        if seed is not None:
            np.random.seed(seed)
        orig_pos = atoms.get_positions()
        orig_cell = atoms.get_cell()
        dpos = np.random.randn(*orig_pos.shape) * pos_shake
        dcell_matrix = np.eye(3)+np.random.randn(3,3)*cell_shake   # Shake matrix to displace volume 
        new_cell = np.dot(orig_cell, dcell_matrix)
        new_pos = orig_pos+dpos
        new_atoms = atoms.copy()
        new_atoms.set_positions(new_pos)
        new_atoms.set_cell(new_cell,scale_atoms=True)
        
        return new_atoms