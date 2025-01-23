# check the counterparts for each glass composition
import json
from pymatgen.ext.matproj import MPRester #, Element
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, GrandPotentialPhaseDiagram, GrandPotPDEntry, CompoundPhaseDiagram
#from pynter import SETTINGS
#from pynter.tools.utils import save_object_as_json, get_object_from_json
import requests
from pymatgen.core.composition import Composition
from pymatgen.io.vasp.outputs import Outcar
from ase.io import read, write
from ase.io.vasp import read_vasp, write_vasp 
from pymatgen.ext.matproj import MPRester
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PatchedPhaseDiagram, CompoundPhaseDiagram, ReactionDiagram
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.reaction_calculator import ComputedReaction

from ase.spacegroup import get_spacegroup
from matplotlib import pyplot as plt
plt.style.use(['vaibhz-sci','no-latex','high-vis','vibrant'])

from ase.io import read
from ase.db import connect
from ase.io.vasp import read_vasp, write_vasp

from ase.io import read
from ase.db import connect
from pymatgen.io.ase import AseAtomsAdaptor
import json
from pymatgen.ext.matproj import MPRester #, Composition, Element
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, GrandPotentialPhaseDiagram, GrandPotPDEntry, CompoundPhaseDiagram
#from pynter import SETTINGS
#from pynter.tools.utils import save_object_as_json, get_object_from_json
import requests
from ase.io import read, write
import re
import numpy as np
from ase.formula import Formula
import itertools
from sympy import symbols, Eq, solve, Rational
from ase.formula import Formula
import collections

# add database 

database = []
for i in range(14):
    # 1 database
    x = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/database/pbe_structure_%s_hdf5/pbe_structure_%s/CONTCAR'%(i,i), format='vasp')
    formula = x.get_chemical_formula(empirical=False)
    atoms = x.get_positions()
    y = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/database/pbe_structure_%s_hdf5/pbe_structure_%s/OUTCAR'%(i,i))
    en = y.final_energy
    #print(formula, en)
    if formula == 'Na12O48P4Si8Zr8':
        entry = ComputedEntry(formula, en, entry_id='Ordered')
        database.append(entry)
    elif formula == 'Na41O144P7Si29Zr24':
        entry = ComputedEntry(formula, en, entry_id='HeNa-1')
        database.append(entry)
    elif formula == 'Na41O143P9Si30Zr20':
        entry = ComputedEntry(formula, en, entry_id='HeNa-2')
        database.append(entry)
    else:
        entry = ComputedEntry(formula, en)
        database.append(entry)
        
database_1 = []
for i in range(13):
    # 1 database
    x = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/database/pbe_structure_%s_hdf5/pbe_structure_%s/CONTCAR'%(i,i), format='vasp')
    formula = x.get_chemical_formula(empirical=False)
    atoms = x.get_positions()
    y = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/database/pbe_structure_%s_hdf5/pbe_structure_%s/OUTCAR'%(i,i))
    en = y.final_energy
    #print(formula, en)
    entry = ComputedEntry(formula, en)
    database_1.append(entry)
    
database_stable = []

for i in range(51):
    # 1 database
    x = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/database_stable/pbe_structure_%s_hdf5/pbe_structure_%s/CONTCAR'%(i,i), format='vasp')
    formula = x.get_chemical_formula(empirical=False)
    atoms = x.get_positions()
    y = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/database_stable/pbe_structure_%s_hdf5/pbe_structure_%s/OUTCAR'%(i,i))
    en = y.final_energy
    #print(formula, en)
    
    if formula == 'Na12O48P4Si8Zr8':
        entry = ComputedEntry(formula, en, entry_id='NaSi-MP-Stable')
        database_stable.append(entry)
        #print(en)
    else:
        entry = ComputedEntry(formula, en)
        database_stable.append(entry)
        
database_glass = []

for i in range(6):
    # 1 database
    x = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/database_glass/pbe_structure_%s_hdf5/pbe_structure_%s/CONTCAR'%(i,i), format='vasp')
    formula = x.get_chemical_formula(empirical=False)
    formula_2 = x.get_chemical_formula(empirical=True)
    atoms = x.get_positions()
    y = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/database_glass/pbe_structure_%s_hdf5/pbe_structure_%s/OUTCAR'%(i,i))
    en = y.final_energy
    #print(formula_2, en)
    if formula_2 == 'Na6O14P2Si3':
        entry = ComputedEntry(formula, en, entry_id='3:2:3')
        database_glass.append(entry)
    elif formula_2 == 'Na6O19P4Si3':
        entry = ComputedEntry(formula, en, entry_id='3:1:3')
        database_glass.append(entry)
    elif formula_2 == 'Na10O19P4Si2':
        entry = ComputedEntry(formula, en, entry_id='5:2:2')
        database_glass.append(entry)
    elif formula_2 == 'Na12O23P2Si6':
        entry = ComputedEntry(formula, en, entry_id='6:1:6')
        database_glass.append(entry)
    elif formula_2 == 'Na12O17P2Si3':
        entry = ComputedEntry(formula, en, entry_id='6:1:3')
        database_glass.append(entry)
    elif formula_2 == 'Na3O8PSi2':
        entry = ComputedEntry(formula, en, entry_id='OQMD-glass')
        database_glass.append(entry)
    else:
        print('No entry')
        
database_nasi = []
for i in range(1):
    # 1 database
    x = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/nasicon_phase_space_new/pbe_structure_%s_hdf5/pbe_structure_%s/CONTCAR'%(i,i), format='vasp')
    formula = x.get_chemical_formula(empirical=False)
    atoms = x.get_positions()
    symmetry = get_spacegroup(x,symprec=1e-5).symbol
    y = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/nasicon_phase_space_new/pbe_structure_%s_hdf5/pbe_structure_%s/OUTCAR'%(i,i))
    en = y.final_energy
    #print(formula, en, symmetry)
    entry = ComputedEntry(formula, en)
    database_nasi.append(entry)

database_nasi_unstable = []
for i in range(18):
    # 1 database
    if i == 2 or i == 8 or i == 11 or i == 17:
        continue
    else:
        x = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/database_unstable/pbe_structure_%s_hdf5/pbe_structure_%s/CONTCAR'%(i,i), format='vasp')
        formula = x.get_chemical_formula(empirical=False)
        atoms = x.get_positions()
        symmetry = get_spacegroup(x,symprec=1e-5).symbol
        y = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/database_unstable/pbe_structure_%s_hdf5/pbe_structure_%s/OUTCAR'%(i,i))
        en = y.final_energy
        #print(formula, en, symmetry)
        if symmetry == 'C c' or symmetry == 'C 2':
            #print(formula, en, symmetry)
            entry = ComputedEntry(formula, en, entry_id='Monoclinic')
            database_nasi_unstable.append(entry)
        else:    
            entry = ComputedEntry(formula, en)
            database_nasi_unstable.append(entry)
            
dt_b_1_5_r = []
for i in range(0,1,1):
    x_1 = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/other_nasicon_compositions_1/nasi_x_1_5_random_hdf5/nasi_x_1_5_random/CONTCAR', format='vasp')
    y_1 = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/other_nasicon_compositions_1/nasi_x_1_5_random_hdf5/nasi_x_1_5_random/OUTCAR')
    formula = x_1.get_chemical_formula(empirical=False)
    atoms = x_1.get_positions()
    en = y_1.final_energy
    datab_p = ComputedEntry(formula, en, entry_id='nasi_x_1_5_ran')
    dt_b_1_5_r.append(datab_p)
    
dt_b_2_p = []
for i in range(0,1,1):
    x_3 = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/other_nasicon_compositions_1/nasi_x_2_planar_hdf5/nasi_x_2_planar/CONTCAR', format='vasp')
    y_3 = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/other_nasicon_compositions_1/nasi_x_2_planar_hdf5/nasi_x_2_planar/OUTCAR')
    formula = x_3.get_chemical_formula(empirical=False)
    atoms = x_3.get_positions()
    en = y_3.final_energy
    datab_2_p = ComputedEntry(formula, en, entry_id='nasi_x_2_planar')
    dt_b_2_p.append(datab_2_p)
    
dt_b_2_r = []
for i in range(0,1,1):
    x_4 = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/other_nasicon_compositions_1/nasi_2_rand_0_hdf5/nasi_2_rand_0/CONTCAR', format='vasp')
    y_4 = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/other_nasicon_compositions_1/nasi_2_rand_0_hdf5/nasi_2_rand_0/OUTCAR')
    formula = x_4.get_chemical_formula(empirical=False)
    atoms = x_4.get_positions()
    en = y_4.final_energy
    datab_2_r = ComputedEntry(formula, en, entry_id='nasi_x_2_rand')
    dt_b_2_r.append(datab_2_r)
    
nasi_ex_0 = []

for i in range(0,1,1):
    x_4 = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_0_hdf5/miss_n_nasi_0/CONTCAR', format='vasp')
    y_4 = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_0_hdf5/miss_n_nasi_0/OUTCAR')
    formula = x_4.get_chemical_formula(empirical=False)
    atoms = x_4.get_positions()
    en = y_4.final_energy
    datab_2_r = ComputedEntry(formula, en, entry_id='nasi_ex_0')
    nasi_ex_0.append(datab_2_r)

nasi_ex_1 = []
for i in range(0,1,1):
    x_4 = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_1_hdf5/miss_n_nasi_1/CONTCAR', format='vasp')
    y_4 = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_1_hdf5/miss_n_nasi_1/OUTCAR')
    formula = x_4.get_chemical_formula(empirical=False)
    atoms = x_4.get_positions()
    en = y_4.final_energy
    datab_2_r = ComputedEntry(formula, en, entry_id='nasi_ex_1')
    nasi_ex_1.append(datab_2_r)
    
nasi_ex_2 = []
for i in range(0,1,1):
    x_4 = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_2_hdf5/miss_n_nasi_2/CONTCAR', format='vasp')
    y_4 = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_2_hdf5/miss_n_nasi_2/OUTCAR')
    formula = x_4.get_chemical_formula(empirical=False)
    atoms = x_4.get_positions()
    en = y_4.final_energy
    datab_2_r = ComputedEntry(formula, en, entry_id='nasi_ex_2')
    nasi_ex_2.append(datab_2_r)
    
nasi_ex_4 = []
for i in range(0,1,1):
    x_4 = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_4_hdf5/miss_n_nasi_4/CONTCAR', format='vasp')
    y_4 = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_4_hdf5/miss_n_nasi_4/OUTCAR')
    formula = x_4.get_chemical_formula(empirical=False)
    atoms = x_4.get_positions()
    en = y_4.final_energy
    datab_2_r = ComputedEntry(formula, en, entry_id='nasi_ex_4')
    nasi_ex_4.append(datab_2_r)
    
nasi_ex_5 = [] 
for i in range(0,1,1):
    x_4 = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_5_hdf5/miss_n_nasi_5/CONTCAR', format='vasp')
    y_4 = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_5_hdf5/miss_n_nasi_5/OUTCAR')
    formula = x_4.get_chemical_formula(empirical=False)
    atoms = x_4.get_positions()
    en = y_4.final_energy
    datab_2_r = ComputedEntry(formula, en, entry_id='nasi_ex_5')
    nasi_ex_5.append(datab_2_r)
    
nasi_ex_7 = []
for i in range(0,1,1):
    x_4 = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_7_hdf5/miss_n_nasi_7/CONTCAR', format='vasp')
    y_4 = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_7_hdf5/miss_n_nasi_7/OUTCAR')
    formula = x_4.get_chemical_formula(empirical=False)
    atoms = x_4.get_positions()
    en = y_4.final_energy
    datab_2_r = ComputedEntry(formula, en, entry_id='nasi_ex_7')
    nasi_ex_7.append(datab_2_r)   
    
nasi_ex_10 = []
for i in range(0,1,1):
    x_4 = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_10_hdf5/miss_n_nasi_10/CONTCAR', format='vasp')
    y_4 = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/extra_nasi_compo/miss_n_nasi_10_hdf5/miss_n_nasi_10/OUTCAR')
    formula = x_4.get_chemical_formula(empirical=False)
    atoms = x_4.get_positions()
    en = y_4.final_energy
    datab_2_r = ComputedEntry(formula, en, entry_id='nasi_ex_10')
    nasi_ex_10.append(datab_2_r)


# Metropolis implementations of the DFT calculations
import os

database_mc = []

base_dir = "/nfshome/deshmukh/vaibhav/NaSICON_dft/database_matropolis_extra_nasi"

# Loop through each subdirectory in the base directory
for dir_name in os.listdir(base_dir):
    dir_path = os.path.join(base_dir, dir_name)
    
    # Check if it's a directory
    if os.path.isdir(dir_path):
        outcar_path = os.path.join(dir_path, "OUTCAR")
        
        # Check if the OUTCAR file exists
        if os.path.isfile(outcar_path):
            print(f"Reading {outcar_path}")
            try:
                x = read(os.path.join(dir_path, "CONTCAR"), format='vasp')
                formula = x.get_chemical_formula(empirical=False)
                atoms = x.get_positions()
                structure = AseAtomsAdaptor.get_structure(x)
                sp_analyzer = SpacegroupAnalyzer(structure)
                sp_grp = sp_analyzer.get_space_group_symbol()
                
                y = Outcar(outcar_path)
                en = y.final_energy
                print(formula, en, sp_grp)
                entry = ComputedEntry(formula, en, data={'spacegroup': sp_grp})
                database_mc.append(entry)                
            except Exception as e:
                print(f"Failed to read {outcar_path}: {e}")
        else:
            print(f"OUTCAR file not found in {dir_path}")
    
    
database_ex = nasi_ex_0 + nasi_ex_1 + nasi_ex_2 + nasi_ex_4 + nasi_ex_5 + nasi_ex_7 + nasi_ex_10

data_my =  dt_b_1_5_r + dt_b_2_p + dt_b_2_r + database_ex

database_all = database + database_stable + database_glass + database_nasi + data_my + database_nasi_unstable + database_mc

pd_db = PhaseDiagram(database_all)

energy = []
formula_all_data = []

# Loop through the stable entries and save the structures to the database
for entry in database_all:
    material_id = entry.entry_id
    form = entry.composition.reduced_formula
    structure = entry.energy_per_atom
    
    #if form in filter_products:
    energy.append(structure)
    formula_all_data.append(form)
    
combined_data = list(zip(formula_all_data, energy))

from collections import defaultdict

formula_energy_dict = defaultdict(list)

# Populate the dictionary with energies
for formula, energy in combined_data:
    formula_energy_dict[formula].append(energy) # Master to search

from collections import defaultdict

# Iterate through the dictionary and replace the lists with the minimum value
for key, values in formula_energy_dict.items():
    if isinstance(values, list) and len(values) > 0:
        formula_energy_dict[key] = min(values)
        
map_key = {
    'Na4Zr2(SiO4)3': 'Na4Zr2Si3O12',
    'NaZr2(PO4)3': 'NaZr2P3O12',
    'Na6Si3(PO7)2': 'Na6Si3P2O14',
    'Zr(PO3)4': 'ZrP4O12',
    'Na5Zr4Si3(PO8)3': 'Na5Zr4Si3P3O24', 
    'Na4Zr6Si(P2O9)4': 'Na4Zr6SiP8O36',
    'Na2Zr2Si(PO6)2': 'Na2Zr2SiP2O12',
    'Na5Zr4Si3(PO8)3': 'Na5Zr4Si3P3O24',
    'Na8Zr6Si5(PO9)4': 'Na8Zr6Si5P4O36',
    'Na10Zr6Si7(PO18)2': 'Na10Zr6Si7P2O36'   
}

formula_energy_dict_new = {map_key.get(k, k): v for k, v in formula_energy_dict.items()}

#(formula_energy_dict_new['Na4Zr2Si3O12'] - formula_energy_dict_new['Na2ZrSi2O7'])*1000

unique_list = list(set(formula_energy_dict_new.keys()))
#print(len(unique_list), len(formula_energy_dict_new.keys()))

possible_products_1 = [i for i in itertools.combinations(unique_list, 1)]
possible_products_2 = [i for i in itertools.combinations(unique_list, 2)]
possible_products_3 = [i for i in itertools.combinations(unique_list, 3)]
possible_products_4 = [i for i in itertools.combinations(unique_list, 4)]
possible_products_5 = [i for i in itertools.combinations(unique_list, 5)]



def reaction_1(interest='NaZr2P3O12', interest_2=Formula('NaZr2P3O12')):
    
    data_1 = collections.defaultdict(list)
    data_2 = collections.defaultdict(list)

    possible_reaction = []
    not_possible_reaction = []
    reaction_energy = []
    reaction_energy_not = []


    # Define the symbols for the coefficients
    coeff_A, coeff_B = symbols('coeff_A coeff_B')

    # Desired composition (Na3Zr2Si2PO12)
    desired = Formula('%s'%interest).count()

    # List of possible products (i[0] and i[1])  # Replace with your list of possible products
    print('#--------------------Possible Reaction and Reaction energy For 2 Products------------------------#')

    for i in possible_products_2:

        # Compositions of the two reactants (i[0] and i[1]) obtained using ASE Formula
        compo_A_formula = Formula(i[0])
        compo_A = compo_A_formula.count()
        compo_B_formula = Formula(i[1])
        compo_B = compo_B_formula.count()

        if all(key in interest_2 for key in compo_A.keys()) and \
            all(key in interest_2 for key in compo_B.keys()): 

            # The energies are in per atom, we need to multiply by the number of atoms in the formula to get the total energy per formula unit
            form_A = sum(list(Formula(i[0]).count().values()))
            form_B = sum(list(Formula(i[1]).count().values()))
            form_int = sum(list(Formula(interest).count().values()))

            en_A = np.min(formula_energy_dict_new.get(i[0]))*form_A
            en_B = np.min(formula_energy_dict_new.get(i[1]))*form_B
            desired_en = np.min(formula_energy_dict_new.get('%s'%interest))*form_int
        


            #print(en_A, en_B, desired_en)
            #print('Energy of ', i[0], 'is ', en_A)
            # Write down the conservation of mass equations for each element
            equations = []
            for element, count in desired.items():
                equation = Eq(count, coeff_A * compo_A.get(element, 0) + coeff_B * compo_B.get(element, 0))
                equations.append(equation)
        
        
            # Solve the system of equations
            solution = solve(equations, (coeff_A, coeff_B), dict=True)  # Use `dict=True` to get a dictionary of solutions
            #print(i)
            # Check if there are any solutions
            if solution:
                # Get the balanced coefficients from the first solution (assuming it's unique) as fractions
                balanced_coeff_A = solution[0].get(coeff_A, 0)
                balanced_coeff_B = solution[0].get(coeff_B, 0)

                #if balanced_coeff_A == balanced_coeff_B:

                #    continue

                if balanced_coeff_A == 1 and balanced_coeff_B == 0:

                    continue
            
                elif balanced_coeff_A == 0 and balanced_coeff_B == 1:
                
                    continue
            
                elif balanced_coeff_A == 0 and balanced_coeff_B == 0:
                
                    continue    
            
                elif balanced_coeff_A < 0 or balanced_coeff_B < 0:
                
                    continue

                #elif balanced_coeff_A + balanced_coeff_B != 1:

                #    continue
            
                else:

                    # Calculate decomposition energy
                    d_energy = ((balanced_coeff_A * en_A + balanced_coeff_B * en_B - desired_en)/form_int) + (-0.0444)    # eV/atom
                    reaction = f"{balanced_coeff_A}*{i[0]} + {balanced_coeff_B}*{i[1]}"

                    if d_energy < 0:
                        data_1['Reactant'].append(interest)
                        data_1['Product'].append(reaction)
                        data_1['DeltaE'].append(round(d_energy,6))
                        possible_reaction.append(reaction)
                        reaction_energy.append(round(d_energy,6))
                        print(f"{'%s'%interest} -> {balanced_coeff_A}*{i[0]} + {balanced_coeff_B}*{i[1]}, deltaE = {round(d_energy, 6)} eV/atom")
                    elif d_energy > 0:
                        data_2['Reactant'].append(interest)
                        data_2['Product'].append(reaction)
                        data_2['DeltaE'].append(round(d_energy,6))
                        not_possible_reaction.append(reaction)
                        reaction_energy_not.append(round(d_energy,6))
                        #print(f"{'%s'%interest} -> {balanced_coeff_A}*{i[0]} + {balanced_coeff_B}*{i[1]}, deltaE = {balanced_coeff_A*en_A + balanced_coeff_B*en_B - desired_en} eV/atom")


                        #print(f"{'%s'%interest} -> {balanced_coeff_A}*{i[0]} + {balanced_coeff_B}*{i[1]}, deltaE = {balanced_coeff_A*en_A + balanced_coeff_B*en_B - desired_en} eV/atom")

            else:
                pass
            #print(f"No solution found for {i[0]} + {i[1]}")
    print('#------------------------------------------------------------------------------------------------#')

    from pandas import DataFrame
    import pandas 
    pandas.set_option('display.max_rows', 100)
    df = DataFrame(data_1, columns=["Reactant", "Product", "DeltaE"])
    df_2 = DataFrame(data_2, columns=["Reactant", "Product", "DeltaE"])
    df.to_csv('possible_reaction_2_%s.csv'%interest, index=True, header=True)
    df_2.to_csv('not_possible_reaction_2_%s.csv'%interest, index=True, header=True)
    
    
    
def reaction_2(interest='NaZr2P3O12', interest_2=Formula('NaZr2P3O12')):
    
    data_1_1 = collections.defaultdict(list)
    data_2_1 = collections.defaultdict(list)

    possible_reaction_3 = []
    not_possible_reaction_3 = []
    reaction_energy_3 = []
    reaction_energy_not_3 = []

    # Define the symbols for the coefficients
    coeff_A, coeff_B, coeff_C = symbols('coeff_A coeff_B coeff_C')

    # Desired composition (Na3Zr2Si2PO12)
    desired = Formula('%s'%interest).count()

    # List of possible products (i[0], i[1], i[2])  # Replace with your list of possible products
    print('#--------------------Possible Reaction and Reaction energy For 3 Products------------------------#')

    for i in possible_products_3:
        # Compositions of the reactants (i[0], i[1], i[2]) obtained using ASE Formula
        compo_A_formula = Formula(i[0])
        compo_A = compo_A_formula.count()
        compo_B_formula = Formula(i[1])
        compo_B = compo_B_formula.count()
        compo_C_formula = Formula(i[2])
        compo_C = compo_C_formula.count()

        # Skip reactions where any pair of reactants have the same composition
        #if compo_A == compo_B or compo_A == compo_C or compo_B == compo_C:
        #    continue

        if all(key in interest_2 for key in compo_A.keys()) and \
            all(key in interest_2 for key in compo_B.keys()) and \
            all(key in interest_2 for key in compo_C.keys()):    
    
            # Retrieve the minimum energy for each reactant
            form_A = sum(list(Formula(i[0]).count().values()))
            form_B = sum(list(Formula(i[1]).count().values()))
            form_C = sum(list(Formula(i[2]).count().values()))
            form_int = sum(list(Formula(interest).count().values()))

            en_A = np.min(formula_energy_dict_new.get(i[0]))*form_A
            en_B = np.min(formula_energy_dict_new.get(i[1]))*form_B
            en_C = np.min(formula_energy_dict_new.get(i[2]))*form_C
            desired_en = np.min(formula_energy_dict_new.get('%s'%interest))*form_int

            # Write down the conservation of mass equations for each element
            equations = []
            for element, count in desired.items():
                equation = Eq(count, coeff_A * compo_A.get(element, 0) + coeff_B * compo_B.get(element, 0) + coeff_C * compo_C.get(element, 0))
                equations.append(equation)

            # Solve the system of equations
            solution = solve(equations, (coeff_A, coeff_B, coeff_C), dict=True)  # Use `dict=True` to get a dictionary of solutions

            #print("Matrix Equation A * x = b:")
            #for equation in equations:
            #    print(equation)
        
            if solution:
                # Get the balanced coefficients from the first solution (assuming it's unique) as fractions
                balanced_coeff_A = solution[0].get(coeff_A, 0)
                balanced_coeff_B = solution[0].get(coeff_B, 0)
                balanced_coeff_C = solution[0].get(coeff_C, 0)

                if balanced_coeff_B == 1 and balanced_coeff_C == 0 and balanced_coeff_A == 0:

                    continue

                elif balanced_coeff_A == 0 or balanced_coeff_B == 1 or balanced_coeff_C == 0:

                    continue

                elif balanced_coeff_A == 0 or balanced_coeff_B == 0 or balanced_coeff_C == 1:

                    continue
        
                elif balanced_coeff_A < 0 or balanced_coeff_B < 0 or balanced_coeff_C < 0:
            
                    continue

                #elif balanced_coeff_A + balanced_coeff_B + balanced_coeff_C != 1:

                #    continue

                else:

                    # Calculate decomposition energy
                    d_energy = ((balanced_coeff_A * en_A + balanced_coeff_B * en_B + balanced_coeff_C * en_C - desired_en)/form_int) + (-0.0444)    # eV/atom
                    reaction = f"{balanced_coeff_A}*{i[0]} + {balanced_coeff_B}*{i[1]} + {balanced_coeff_C}*{i[2]}"
            
                    if d_energy < 0:
                        data_1_1['Reactant'].append(interest)
                        data_1_1['Product'].append(reaction)
                        data_1_1['DeltaE'].append(round(d_energy, 20))
                        possible_reaction_3.append(reaction)
                        reaction_energy_3.append(round(d_energy, 20))
                        print(f"{'%s'%interest} -> {balanced_coeff_A}*{i[0]} + {balanced_coeff_B}*{i[1]} + {balanced_coeff_C}*{i[2]}, deltaE = {round(d_energy, 20)} eV/atom")
                    elif d_energy > 0:
                        data_2_1['Reactant'].append(interest)
                        data_2_1['Product'].append(reaction)
                        data_2_1['DeltaE'].append(round(d_energy, 20))
                        not_possible_reaction_3.append(reaction)
                        reaction_energy_not_3.append(round(d_energy, 20))
        else:
            pass

    print('#------------------------------------------------------------------------------------------------#')

    from pandas import DataFrame
    import pandas 
    pandas.set_option('display.max_rows', 100)
    df_1 = DataFrame(data_1_1, columns=["Reactant", "Product", "DeltaE"])
    df_2_1 = DataFrame(data_2_1, columns=["Reactant", "Product", "DeltaE"])
    df_1.to_csv('possible_reaction_3_%s.csv'%interest, index=True, header=True)
    df_2_1.to_csv('not_possible_reaction_3_%s.csv'%interest, index=True, header=True)
    

def reaction_3(interest='NaZr2P3O12', interest_2=Formula('NaZr2P3O12')):
    
    data_1_2 = collections.defaultdict(list)
    data_2_2 = collections.defaultdict(list)


    possible_reaction_4 = []
    not_possible_reaction_4 = []
    reaction_energy_4 = []
    reaction_energy_not_4 = []

    # Define the symbols for the coefficients
    coeff_A, coeff_B, coeff_C, coeff_D = symbols('coeff_A coeff_B coeff_C coeff_D')

    # Desired composition (Na3Zr2Si2PO12)
    desired = Formula('%s'%interest).count()

    # List of possible products (i[0], i[1], i[2])  # Replace with your list of possible products
    print('#--------------------Possible Reaction and Reaction energy For 4 Products------------------------#')

    for i in possible_products_4:
        # Compositions of the reactants (i[0], i[1], i[2]) obtained using ASE Formula
        compo_A_formula = Formula(i[0])
        compo_A = compo_A_formula.count()
        compo_B_formula = Formula(i[1])
        compo_B = compo_B_formula.count()
        compo_C_formula = Formula(i[2])
        compo_C = compo_C_formula.count()
        compo_D_formula = Formula(i[3])
        compo_D = compo_D_formula.count()

        if all(key in interest_2 for key in compo_A.keys()) and \
            all(key in interest_2 for key in compo_B.keys()) and \
            all(key in interest_2 for key in compo_C.keys()) and \
            all(key in interest_2 for key in compo_D.keys()):

            # Retrieve the minimum energy for each reactant
            form_A = sum(list(Formula(i[0]).count().values()))
            form_B = sum(list(Formula(i[1]).count().values()))
            form_C = sum(list(Formula(i[2]).count().values()))
            form_D = sum(list(Formula(i[3]).count().values()))
            form_int = sum(list(Formula(interest).count().values()))

            en_A = np.min(formula_energy_dict_new.get(i[0]))*form_A
            en_B = np.min(formula_energy_dict_new.get(i[1]))*form_B
            en_C = np.min(formula_energy_dict_new.get(i[2]))*form_C
            en_D = np.min(formula_energy_dict_new.get(i[3]))*form_D

            desired_en = np.min(formula_energy_dict_new.get('%s'%interest))*form_int

            # Write down the conservation of mass equations for each element
            equations = []
            for element, count in desired.items():
                equation = Eq(count, coeff_A * compo_A.get(element, 0) + coeff_B * compo_B.get(element, 0) + coeff_C * compo_C.get(element, 0) + coeff_D * compo_D.get(element, 0))
                equations.append(equation)

            # Solve the system of equations
            solution = solve(equations, (coeff_A, coeff_B, coeff_C, coeff_D), dict=True)  # Use `dict=True` to get a dictionary of solutions

            if solution:
                # Get the balanced coefficients from the first solution (assuming it's unique) as fractions
                balanced_coeff_A = solution[0].get(coeff_A, 0)
                balanced_coeff_B = solution[0].get(coeff_B, 0)
                balanced_coeff_C = solution[0].get(coeff_C, 0)
                balanced_coeff_D = solution[0].get(coeff_D, 0)

                if balanced_coeff_B == 1 and balanced_coeff_C == 0 and balanced_coeff_A == 0 and balanced_coeff_D == 0:

                    continue

                elif balanced_coeff_A == 0 or balanced_coeff_B == 1 or balanced_coeff_C == 0 or balanced_coeff_D == 0:

                    continue

                elif balanced_coeff_A == 0 or balanced_coeff_B == 0 or balanced_coeff_C == 1 or balanced_coeff_D == 0:

                    continue

                elif balanced_coeff_A == 0 or balanced_coeff_B == 0 or balanced_coeff_C == 0 or balanced_coeff_D == 1:

                    continue
        
                elif balanced_coeff_A < 0 or balanced_coeff_B < 0 or balanced_coeff_C < 0 or balanced_coeff_D < 0:
            
                    continue

                else:

                    # Calculate decomposition energy
                    d_energy = ((balanced_coeff_A * en_A + balanced_coeff_B * en_B + balanced_coeff_C * en_C + balanced_coeff_D * en_D - desired_en)/form_int) + (-0.0444)  # eV/atom
                    reaction = f"{balanced_coeff_A}*{i[0]} + {balanced_coeff_B}*{i[1]} + {balanced_coeff_C}*{i[2]} + {balanced_coeff_D}*{i[3]}"
 
                    if d_energy < 0:
                        data_1_2['Reactant'].append(interest)
                        data_1_2['Product'].append(reaction)
                        data_1_2['DeltaE'].append(round(d_energy, 20))
                        possible_reaction_4.append(reaction)
                        reaction_energy_4.append(round(d_energy, 20))
                        print(f"{'%s'%interest} -> {balanced_coeff_A}*{i[0]} + {balanced_coeff_B}*{i[1]} + {balanced_coeff_C}*{i[2]} + {balanced_coeff_D}*{i[3]}, deltaE = {round(d_energy, 20)} eV/atom")
                    elif d_energy > 0:
                        data_2_2['Reactant'].append(interest)
                        data_2_2['Product'].append(reaction)
                        data_2_2['DeltaE'].append(round(d_energy, 20))
                        not_possible_reaction_4.append(reaction)
                        reaction_energy_not_4.append(round(d_energy, 20))
        else:
            pass

    print('#------------------------------------------------------------------------------------------------#')

    from pandas import DataFrame
    import pandas 
    pandas.set_option('display.max_rows', 100)
    
    df_1_3 = DataFrame(data_1_2, columns=["Reactant", "Product", "DeltaE"])
    df_2_3 = DataFrame(data_2_2, columns=["Reactant", "Product", "DeltaE"])
    
    df_1_3.to_csv('possible_reaction_4_%s.csv'%interest, index=True, header=True)
    df_2_3.to_csv('not_possible_reaction_4_%s.csv'%interest, index=True, header=True)
    
def reaction_4(interest='NaZr2P3O12', interest_2=Formula('NaZr2P3O12')):
    
    data_1_4 = collections.defaultdict(list)
    data_2_4 = collections.defaultdict(list)

    possible_reaction_5 = []
    not_possible_reaction_5 = []
    reaction_energy_5 = []
    reaction_energy_not_5 = []

    # Define the symbols for the coefficients
    coeff_A, coeff_B, coeff_C, coeff_D, coeff_E = symbols('coeff_A coeff_B coeff_C coeff_D coeff_E')

    # Desired composition (Na3Zr2Si2PO12)
    desired = Formula('%s'%interest).count()

    # List of possible products (i[0], i[1], i[2])  # Replace with your list of possible products
    print('#--------------------Possible Reaction and Reaction energy For 5 Products------------------------#')

    for i in possible_products_5:
        # Compositions of the reactants (i[0], i[1], i[2]) obtained using ASE Formula
        compo_A_formula = Formula(i[0])
        compo_A = compo_A_formula.count()
        compo_B_formula = Formula(i[1])
        compo_B = compo_B_formula.count()
        compo_C_formula = Formula(i[2])
        compo_C = compo_C_formula.count()
        compo_D_formula = Formula(i[3])
        compo_D = compo_D_formula.count()
        compo_E_formula = Formula(i[4])
        compo_E = compo_E_formula.count()

        if all(key in interest_2 for key in compo_A.keys()) and \
            all(key in interest_2 for key in compo_B.keys()) and \
            all(key in interest_2 for key in compo_C.keys()) and \
            all(key in interest_2 for key in compo_D.keys()) and \
            all(key in interest_2 for key in compo_E.keys()):

            # Retrieve the minimum energy for each reactant
            form_A = sum(list(Formula(i[0]).count().values()))
            form_B = sum(list(Formula(i[1]).count().values()))
            form_C = sum(list(Formula(i[2]).count().values()))
            form_D = sum(list(Formula(i[3]).count().values()))
            form_E = sum(list(Formula(i[4]).count().values()))
            form_int = sum(list(Formula(interest).count().values()))

            en_A = np.min(formula_energy_dict_new.get(i[0]))*form_A
            en_B = np.min(formula_energy_dict_new.get(i[1]))*form_B
            en_C = np.min(formula_energy_dict_new.get(i[2]))*form_C
            en_D = np.min(formula_energy_dict_new.get(i[3]))*form_D
            en_E = np.min(formula_energy_dict_new.get(i[4]))*form_E

            desired_en = np.min(formula_energy_dict_new.get('%s'%interest))*form_int

            # Write down the conservation of mass equations for each element
            equations = []
            for element, count in desired.items():
                equation = Eq(count, coeff_A * compo_A.get(element, 0) + coeff_B * compo_B.get(element, 0) + coeff_C * compo_C.get(element, 0) + coeff_D * compo_D.get(element, 0) + coeff_E * compo_E.get(element, 0))
                equations.append(equation)

            # Solve the system of equations
            solution = solve(equations, (coeff_A, coeff_B, coeff_C, coeff_D, coeff_E), dict=True)  # Use `dict=True` to get a dictionary of solutions

            if solution:
                # Get the balanced coefficients from the first solution (assuming it's unique) as fractions
                balanced_coeff_A = solution[0].get(coeff_A, 0)
                balanced_coeff_B = solution[0].get(coeff_B, 0)
                balanced_coeff_C = solution[0].get(coeff_C, 0)
                balanced_coeff_D = solution[0].get(coeff_D, 0)
                balanced_coeff_E = solution[0].get(coeff_E, 0)
    

                if balanced_coeff_B == 1 and balanced_coeff_C == 0 and balanced_coeff_A == 0 and balanced_coeff_D == 0 and balanced_coeff_E == 0:

                    continue

                elif balanced_coeff_A == 0 or balanced_coeff_B == 1 or balanced_coeff_C == 0 or balanced_coeff_D == 0 or balanced_coeff_E == 0:

                    continue

                elif balanced_coeff_A == 0 or balanced_coeff_B == 0 or balanced_coeff_C == 1 or balanced_coeff_D == 0 or balanced_coeff_E == 0:

                    continue

                elif balanced_coeff_A == 0 or balanced_coeff_B == 0 or balanced_coeff_C == 0 or balanced_coeff_D == 1 or balanced_coeff_E == 0:

                    continue
        
                elif balanced_coeff_A == 0 or balanced_coeff_B == 0 or balanced_coeff_C == 0 or balanced_coeff_D == 0 or balanced_coeff_E == 1:

                    continue
            
                elif balanced_coeff_A < 0 or balanced_coeff_B < 0 or balanced_coeff_C < 0 or balanced_coeff_D < 0 or balanced_coeff_E < 0:
            
                    continue

                else:

                    # Calculate decomposition energy
                    d_energy = ((balanced_coeff_A * en_A + balanced_coeff_B * en_B + balanced_coeff_C * en_C + balanced_coeff_D * en_D + balanced_coeff_E * en_E - desired_en)/form_int) + (-0.0444)    # eV/atom
                    reaction = f"{interest} -> {balanced_coeff_A}*{i[0]} + {balanced_coeff_B}*{i[1]} + {balanced_coeff_C}*{i[2]} + {balanced_coeff_D}*{i[3]} + {balanced_coeff_E}*{i[4]}"
 
                    if d_energy < 0:
                        data_1_4['Reactant'].append(interest)
                        data_1_4['Product'].append(reaction)
                        data_1_4['DeltaE'].append(round(d_energy, 3))
                        possible_reaction_5.append(reaction)
                        reaction_energy_5.append(round(d_energy, 3))
                        print(f"{'%s'%interest} -> {balanced_coeff_A}*{i[0]} + {balanced_coeff_B}*{i[1]} + {balanced_coeff_C}*{i[2]} + {balanced_coeff_D}*{i[3]} + {balanced_coeff_E}*{i[4]}, deltaE = {round(d_energy, 3)} eV/atom")
                    elif d_energy > 0:
                        data_2_4['Reactant'].append(interest)
                        data_2_4['Product'].append(reaction)
                        data_2_4['DeltaE'].append(round(d_energy, 3))
                        not_possible_reaction_5.append(reaction)
                        reaction_energy_not_5.append(round(d_energy, 3))
        else:
            pass

    print('#------------------------------------------------------------------------------------------------#')

    from pandas import DataFrame
    import pandas 
    pandas.set_option('display.max_rows', 100)
    df_1_4 = DataFrame(data_1_4, columns=["Reactant", "Product", "DeltaE"])
    df_2_4 = DataFrame(data_2_4, columns=["Reactant", "Product", "DeltaE"])
    df_1_4.to_csv('possible_reaction_5_%s.csv'%interest, index=True, header=True)
    df_2_4.to_csv('not_possible_reaction_5_%s.csv'%interest, index=True, header=True)
    
#list_nasi = ['NaZr2P3O12', 'Na5Zr6Si2P7O36', 'Na3Zr2Si2PO12', 'Na5Zr4Si3P3O24','Na41Zr24Si29P7O144',
#             'Na7Zr4Si5PO24', 'Na2Zr2SiP2O12', 'Na10Zr6Si7P2O36', 'Na8Zr6Si5P4O36', 
#             'Na3Zr4SiP5O24', 'Na11Zr8Si7P5O48', 'Na4Zr6SiP8O36', 'Na13Zr12Si7P11O72', 
#             'Na7Zr12SiP17O72', 'Na4Zr2Si3O12']
#list_nasi = ['Na11Zr8Si7P5O48', 'Na4Zr6SiP8O36', 'Na13Zr12Si7P11O72', 'Na7Zr12SiP17O72', 'Na4Zr2Si3O12']
#list_nasi = ['Na2ZrSi2O7', 'Na2ZrSiO5', 'Na2ZrSi4O11']  
list_nasi = ['Na41Zr20Si39P9O143']

for i in range(len(list_nasi)):
    reaction_1(interest=list_nasi[i], interest_2=Formula(list_nasi[i])) 
    reaction_2(interest=list_nasi[i], interest_2=Formula(list_nasi[i]))
    reaction_3(interest=list_nasi[i], interest_2=Formula(list_nasi[i]))
    reaction_4(interest=list_nasi[i], interest_2=Formula(list_nasi[i]))

#reaction_1(interest=list_nasi[0], interest_2=Formula(list_nasi[0]))
#reaction_2(interest=list_nasi[0], interest_2=Formula(list_nasi[0]))
#reaction_3(interest=list_nasi[0], interest_2=Formula(list_nasi[0]))
#reaction_4(interest=list_nasi[0], interest_2=Formula(list_nasi[0]))
