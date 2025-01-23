# check the counterparts for each glass composition
import json
from pymatgen.ext.matproj import MPRester, Element
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, GrandPotentialPhaseDiagram, GrandPotPDEntry, CompoundPhaseDiagram
#from pynter import SETTINGS
#from pynter.tools.utils import save_object_as_json, get_object_from_json
import requests
from pymatgen.core.composition import Composition
from pymatgen.io.vasp.outputs import Outcar
from ase.io import read, write
from ase.io.vasp import read_vasp, write_vasp 
from pymatgen.ext.matproj import MPRester, Composition 
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
from pymatgen.ext.matproj import MPRester, Composition, Element
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



##############################

interest = 'NaZr2P3O12'

interest_2 = Formula(interest)

#############################


database = []
for i in range(14):
    # 1 database
    x = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/database/pbe_structure_%s_hdf5/pbe_structure_%s/CONTCAR'%(i,i), format='vasp')
    formula = x.get_chemical_formula(empirical=False)
    atoms = x.get_positions()
    y = Outcar('/nfshome/deshmukh/vaibhav/NaSICON_dft/database/pbe_structure_%s_hdf5/pbe_structure_%s/OUTCAR'%(i,i))
    en = y.final_energy
    print(formula, en)
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
    print(formula, en)
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
        print(en)
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
    print(formula, en, symmetry)
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

data_my =  dt_b_1_5_r + dt_b_2_p + dt_b_2_r

database_all = database + database_stable + database_glass + database_nasi + data_my + database_nasi_unstable 

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
}

formula_energy_dict_new = {map_key.get(k, k): v for k, v in formula_energy_dict.items()}

unique_list = list(set(formula_energy_dict_new.keys()))
print(len(unique_list), len(formula_energy_dict_new.keys()))
unique_list

possible_products_1 = [i for i in itertools.combinations(unique_list, 1)]
possible_products_2 = [i for i in itertools.combinations(unique_list, 2)]
possible_products_3 = [i for i in itertools.combinations(unique_list, 3)]
possible_products_4 = [i for i in itertools.combinations(unique_list, 4)]
possible_products_5 = [i for i in itertools.combinations(unique_list, 5)]

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
                d_energy = (balanced_coeff_A * en_A + balanced_coeff_B * en_B + balanced_coeff_C * en_C + balanced_coeff_D * en_D + balanced_coeff_E * en_E - desired_en)/form_int    # eV/atom
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
