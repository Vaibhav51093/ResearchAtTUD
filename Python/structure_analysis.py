# The script to analyse local structure 
# ASE python follwoing packages 
# 1. Neightbourlist 
# 2. Connectivity matrix  
# 3. Bond analysis 
from ase import neighborlist
from ase.data import covalent_radii 
from ase.calculators.neighborlist import NeighborList
from ase.neighborlist import neighbor_list
from ase.io import read, write 
from ase import Atoms 
from scipy import sparse
import numpy as np 
from ase.io.pov import get_bondpairs 

def get_symbols(structure="300k_amorphous_ordered.cif", form="cif"):
    # Import the structure for analysis
    # Get the symbols
    symbols = []  

    strcut = read(structure, format=form) 

    for i in range(len(strcut)):
        symbols.append(strcut.symbols[i])

# Get neighbour list 
def nn(structure, cut):
    cutOff  = cut * covalent_radii[structure.numbers]    # 1.6 initial
    neighborList = neighborlist.NeighborList(cutoffs=cutOff, self_interaction=False, bothways=True)
    neighborList.update(structure)
    return neighborList

# Get index for each elements  
    Zr_index = [strcu.index for strcu in ase_atoms_0 if strcu.symbol == "Zr"]
    O_index = [strcu.index for strcu in ase_atoms_0 if strcu.symbol == "O"]
    Si_index = [strcu.index for strcu in ase_atoms_0 if strcu.symbol == "Si"]
    P_index = [strcu.index for strcu in ase_atoms_0 if strcu.symbol == "P"]
    Na_index = [strcu.index for strcu in ase_atoms_0 if strcu.symbol == "Na"]
    
    nl = nn(structure=ase_atoms_0, cut=3.5)

    indi_na = []
    indi_zr = []
    indi_p = []
    indi_si = []
    indi_o = []
    for i,j,k,l,m in zip(Na_index,Zr_index,O_index,Si_index,P_index):
        indices_na, offsets_na = nl.get_neighbors(i)
        indices_Zr, offsets_Zr = nl.get_neighbors(j)
        indices_O, offsets_O = nl.get_neighbors(k)
        indices_Si, offsets_Si = nl.get_neighbors(l)
        indices_P, offsets_P = nl.get_neighbors(m)
        indi_na.append(indices_na)
        indi_zr.append(indices_Zr)
        indi_o.append(indices_O)
        indi_si.append(indices_Si)
        indi_p.append(indices_P)

    for i, k in zip(Na_index, indi_na):
        for j in k: 
            if ase_atoms_0.symbols[j]=='Na' or ase_atoms_0.symbols[j]=='O':
                dist = ase_atoms_0.get_distance(i,j) 
                #na_o.append(dist)
                if 2.5<dist<2.8:
                    na_o.append(dist)
                elif 2.8<dist<3.20:
                    na_mid_o.append(dist) 

    for i, k in zip(Zr_index, indi_zr):
        for j in k: 
            if ase_atoms_0.symbols[j]=='Zr' or ase_atoms_0.symbols[j]=='O':
                dist = ase_atoms_0.get_distance(i,j) 
                #zr_o.append(dist)
                if 1.97<dist<2.1:
                    zr_o.append(dist)

    for i, k in zip(Si_index, indi_si):
        for j in k: 
            if ase_atoms_0.symbols[j]=='Si' or ase_atoms_0.symbols[j]=='O':
                dist = ase_atoms_0.get_distance(i,j) 
                #si_o.append(dist)
                if 1.50<dist<1.7:
                    si_o.append(dist)

    for i, k in zip(P_index, indi_p):
        for j in k: 
            if ase_atoms_0.symbols[j]=='P' or ase_atoms_0.symbols[j]=='O':
                dist = ase_atoms_0.get_distance(i,j) 
                #p_o.append(dist)
                if 1.50<dist<1.7:
                    p_o.append(dist)

    for i, k in zip(Na_index, indi_na):
        for j in k: 
            if ase_atoms_0.symbols[j]=='Na' or ase_atoms_0.symbols[j]=='Na':
                dist = ase_atoms_0.get_distance(i,j) 
                #na_na.append(dist)
                if dist<3.5:
                    na_na.append(dist)