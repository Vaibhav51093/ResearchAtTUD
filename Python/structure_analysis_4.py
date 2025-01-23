# check the counterparts for each glass composition
from tempfile import tempdir
from pyiron import Project
import numpy as np
import pandas
from jinja2 import Template
import matplotlib.pyplot as plt 
import scipy.constants as sc
from scipy.integrate import cumtrapz
import os
from pyiron import ase_to_pyiron, pyiron_to_ase
from ase.io import read, write
import random
import numpy as np
from ase.atoms import Atoms, Atom
from itertools import product
from ase.formula import Formula
from ase.formula import Formula
from fileinput import filename
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
from ovito.io.ase import ovito_to_ase, ase_to_ovito
import numpy as np 
plt.style.use(['vaibhz-sci','no-latex','high-vis','vibrant'])
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/frank/zr_vacancy_2")               #  
pr_2 = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/frank/zr_vacancy")               #  
pr_0 = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/frank/potetial_phase_diagram")   #  
pr_3 = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/nasi_2_4_random")     #  

def structure_analysis(job_1, job_2):
    
    # Structures to analyze
    pristine_structure = pyiron_to_ase(job_1.get_structure(iteration_step=-1))
    intial_structure = pyiron_to_ase(job_2.get_structure(iteration_step=0))
    final_structure = pyiron_to_ase(job_2.get_structure(iteration_step=-1))
    
    # Indices of atoms in the structures
    na_index_ref = [Atom.index for Atom in pristine_structure if Atom.symbol == 'Na']
    zr_index_ref = [Atom.index for Atom in pristine_structure if Atom.symbol == 'Zr']
    p_index_ref = [Atom.index for Atom in pristine_structure if Atom.symbol == 'P']
    o_index_ref = [Atom.index for Atom in pristine_structure if Atom.symbol == 'O']
    si_index_ref = [Atom.index for Atom in pristine_structure if Atom.symbol == 'Si']

    na_index_init = [Atom.index for Atom in intial_structure if Atom.symbol == 'Na']
    zr_index_init = [Atom.index for Atom in intial_structure if Atom.symbol == 'Zr']
    p_index_init = [Atom.index for Atom in intial_structure if Atom.symbol == 'P']
    o_index_init = [Atom.index for Atom in intial_structure if Atom.symbol == 'O']
    si_index_init = [Atom.index for Atom in intial_structure if Atom.symbol == 'Si']

    na_index_fin = [Atom.index for Atom in final_structure if Atom.symbol == 'Na']
    zr_index_fin = [Atom.index for Atom in final_structure if Atom.symbol == 'Zr']
    p_index_fin = [Atom.index for Atom in final_structure if Atom.symbol == 'P']
    o_index_fin = [Atom.index for Atom in final_structure if Atom.symbol == 'O']
    si_index_fin = [Atom.index for Atom in final_structure if Atom.symbol == 'Si']
    
    print("Total no of Zr atoms in the pristine structure: ", len(zr_index_ref))
    print("Total no of Zr atoms in the initial structure: ", len(zr_index_init))
    print("Total Zr vacancies: ", len(zr_index_ref) - len(zr_index_init))
    
    # Find the missing Zr index in the intial structure, those we reomved in ref structure to create vacancy

    zr_index_ref_n = []
    zr_index_inti_n = []


    for i in zr_index_ref:
    
        pos_1 = pristine_structure.positions[i]
    
        for j in zr_index_init:
        
            pos_2 = intial_structure.positions[j]
        
            if np.allclose(0, pos_1[0]-pos_2[0], atol=0.0001): # If two index are close to each other, then append the index

                #print(pos_1[0]-pos_2[0])
                zr_index_ref_n.append(i)
                zr_index_inti_n.append(j)
            
            else:
            
                continue
        
    # Unique index of Zr in prestine and intial structure

    zr_index_ref_n = list(set(zr_index_ref_n))
    zr_index_inti_n = list(set(zr_index_inti_n))

    zr_index_ref_vac  = [i for i in zr_index_ref if i not in zr_index_ref_n]


    # Find the index of Na in the intial structure, those we added in ref structure to compenstate the vacancy

    na_index_ref_n = []
    na_index_inti_n = []

    for i in na_index_ref:
    
        pos_1 = pristine_structure.positions[i]
    
        for j in na_index_init:
        
            pos_2 = intial_structure.positions[j]
        
            if np.allclose(0, pos_1[0]-pos_2[0], atol=0.0001): # If two index are close to each other, then append the index

                #print(pos_1[0]-pos_2[0])
                na_index_ref_n.append(i)
                na_index_inti_n.append(j)
            
            else:
            
                continue
        
    # Unique index of Na in prestine and intial structure

    na_index_ref_n = list(set(na_index_ref_n))
    na_index_inti_n = list(set(na_index_inti_n))

    na_index_inti_vac  = [i for i in na_index_init if i not in na_index_inti_n]

    # Find the index of O in the intial structure, those we reomved in ref structure to create vacancy

    o_index_ref_n = []
    o_index_inti_n = []

    for i in o_index_ref:
    
        pos_1 = pristine_structure.positions[i]
    
        for j in o_index_init:
        
            pos_2 = intial_structure.positions[j]
        
            if np.allclose(0, pos_1[0]-pos_2[0], atol=0.0001): # If two index are close to each other, then append the index

                #print(pos_1[0]-pos_2[0])
                o_index_ref_n.append(i)
                o_index_inti_n.append(j)
            
            else:
            
                continue
        
    # Unique index of O in prestine and intial structure

    o_index_ref_n = list(set(o_index_ref_n))
    o_index_inti_n = list(set(o_index_inti_n))

    o_index_ref_vac  = [i for i in o_index_ref if i not in o_index_ref_n]

    print('Total Removed Zr atoms: ', len(zr_index_ref_vac))
    print('Total Added Na atoms: ', len(na_index_inti_vac))
    print('Total Removed O atoms: ', len(o_index_ref_vac))
    
    from ase import neighborlist as nl
    from ase import Atoms
    from scipy import sparse
    
    from ase import Atoms

    class NeighborList:
    
        def __init__(self, atoms, cutoff):
            self.atoms = atoms
            self.cutoff = cutoff
            self.neighbor_list = {}
            self.neighbor_list_symbol = {}

            # Build the neighbor list
            for i in range(len(self.atoms)):
                neighbors = []
                symbol = []
                for j in range(len(self.atoms)):
                    if self.atoms.get_distance(i, j) < self.cutoff:  # Pass atom indices
                        neighbors.append(j)
                        symbol.append(self.atoms[j].symbol)
                self.neighbor_list[i] = neighbors
                self.neighbor_list_symbol[i] = symbol

        def get_neighbors(self, index):
        
            # remove the index of the atom itself
            self.neighbor_list[index].remove(index)
            self.neighbor_list_symbol[index].remove(self.atoms[index].symbol)
            return self.neighbor_list[index], self.neighbor_list_symbol[index]
        
    cn_na_o_init = []          # Coodination number of Na and O where Na is at vacancy site
    cn_na_o_final = []         # Coodination number of Na and O where Na is at vacancy site
    cn_na_o_index = []         # Index of Na, where Na is at vacancy site
    cn_o_index_init = []            # Index of O, where Na is at vacancy site and they are surrounded by these O atoms
    cn_o_index_final = []           # Index of O, where Na is at vacancy site and they are surrounded by these O atoms  
        
    for i in na_index_inti_vac:
    
        index, symbol = NeighborList(intial_structure, 3).get_neighbors(i)
        index_2, symbol_2 = NeighborList(final_structure, 3).get_neighbors(i)
    
        o_count = symbol.count('O')
        si_count = symbol.count('Si')
        p_count = symbol.count('P')
        na_count = symbol.count('Na')
    
        o_count_2 = symbol_2.count('O')
        si_count_2 = symbol_2.count('Si')
        p_count_2 = symbol_2.count('P')
        na_count_2 = symbol_2.count('Na')
    
        # Check index intial and final structure

        o_index = []
        o_index_2 = []
    
        for k, l in zip(index, index_2):
        
            if intial_structure[k].symbol == 'O':
            
                o_index.append(k)

            if final_structure[l].symbol == 'O':
                
                o_index_2.append(l)
    
        common_elements = []

        for item in o_index:
            if item in o_index_2:
                common_elements.append(item)

        count_common = len(common_elements) # No of same O atoms in intial and final structure and also Coodination number
            
        if count_common >= 4:                                                        
            cn_na_o_init.append(len(o_index))
            cn_na_o_final.append(count_common)
            cn_na_o_index.append(i)
            cn_o_index_init.append(o_index)
            cn_o_index_final.append(common_elements)       
                
    return cn_na_o_init,cn_na_o_final, cn_na_o_index, cn_o_index_init, cn_o_index_final    

c_no_init = []
c_no_final = []
c_no_na_index = []
c_no_o_index_init = []
c_no_o_index_final = []

for i in range(1,30,1):
    cn_init, cn_final, na_index, o_index_init, o_index_fin = structure_analysis(pr['nasi_2_5_0'],pr['nasi_2_5_zr_9_%s_4'%i])
    c_no_init.append(cn_init)
    c_no_final.append(cn_final)
    c_no_na_index.append(na_index)
    c_no_o_index_init.append(o_index_init)
    c_no_o_index_final.append(o_index_fin)
    
x = []
x_std = []
y = []
y_std = []
for i, j in zip(c_no_init, c_no_final):
    x.append(np.average(i))
    x_std.append(np.std(i))
    y.append(np.average(j))
    y_std.append(np.std(j))
    
z = range(1,len(x)+1)

data = []

for i in z:
    data.append(i/54)
    
plt.errorbar(data, x,yerr=x_std, fmt='o--', label='Initial')
plt.errorbar(data, y, yerr=y_std, fmt='o--', label='Optimized')
plt.xlabel('y in $\mathrm{Na_{1+x+4y}Zr_{2-y}P_{3-x}Si_{x}O_{12}}$ ')
plt.ylabel('Avg. Coodination number Na-O')
plt.grid(True, axis='y', linestyle='--')
plt.grid(True, axis='x', linestyle='--')
plt.legend()
plt.ylim(2,8)
plt.savefig('cn_na_o_2.png')
plt.show()