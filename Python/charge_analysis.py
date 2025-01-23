#!/usr/bin/env python

import numpy as np
from ase.io import read, write
from ase.units import Bohr
from ase.visualize import view

print('------------------------------------------------------------------------------')

print('Welcome to the charge analysis script using Bader Charges!')
print('This is a official script belongs to Msc. Vaibhav Arun Deshmukh.')

print('------------------------------------------------------------------------------')

def attach_charges(atoms, fileobj='ACF.dat', displacement=1e-4, use_diff=True,
                   use_bohr=True):
    """Attach the charges from the fileobj to the Atoms."""
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    sep = '---------------'
    i = 0 # Counter for the lines
    k = 0 # Counter of sep
    assume6columns = False
    for line in fileobj:
        if line[0] == '\n': # check if there is an empty line in the 
            i -= 1          # head of ACF.dat file
        if i == 0:
            headings = line
            if 'BADER' in headings.split():
                j = headings.split().index('BADER')
            elif 'CHARGE' in headings.split():
                j = headings.split().index('CHARGE')
            else:
                print('Can\'t find keyword "BADER" or "CHARGE".' \
                +' Assuming the ACF.dat file has 6 columns.')
                j = 4
                assume6columns = True
        if sep in line: # Stop at last seperator line
            if k == 1:
                break
            k += 1
        if not i > 1:
            pass
        else:
            words = line.split()
            if assume6columns is True:
                if len(words) != 6:
                    raise IOError('Number of columns in ACF file incorrect!\n'
                                  'Check that Bader program version >= 0.25')
                
            atom = atoms[int(words[0]) - 1]
            if use_diff:
                atom.charge = atom.number - float(words[j])
            else:
                atom.charge = float(words[j])
            if displacement is not None: # check if the atom positions match
                if use_bohr:
                    xyz = np.array([float(w) for w in words[1:4]]) * Bohr
                else:
                    xyz = np.array([float(w) for w in words[1:4]])
                assert np.linalg.norm(atom.position - xyz) < displacement
        i += 1


print('Reading the POSCAR file...')
atoms = read('POSCAR')
print('Reading the ACF.dat file...')



attach_charges(atoms, 'ACF.dat', use_bohr=False, use_diff=False)

# print exchange the charges of individual atoms

for atom in atoms:
    if atom.symbol == 'Na':
        atom.charge = 7 - atom.charge
    elif atom.symbol == 'O':
        atom.charge = 6 - atom.charge
    elif atom.symbol == 'P':
        atom.charge = 5 - atom.charge
    elif atom.symbol == 'Si':
        atom.charge = 4 - atom.charge
    elif atom.symbol == 'Zr':
        atom.charge = 12 - atom.charge

na_c = []
zr_c = []
o_c = []
p_c = []
si_c = []


for i in atoms:
    if i.symbol == 'Na':
        na_c.append(np.round(i.charge, 2))
    elif i.symbol == 'Zr':
        zr_c.append(np.round(i.charge, 2))
    elif i.symbol == 'O':
        o_c.append(np.round(i.charge, 2))
    elif i.symbol == 'P':
        p_c.append(np.round(i.charge, 2))
    elif i.symbol == 'Si':
        si_c.append(np.round(i.charge, 2))

total_c = sum(na_c) + sum(zr_c) + sum(o_c) + sum(p_c) + sum(si_c)
print('----------------------------Charge Analysis (Average)-------------------------')
if sum(na_c) != 0:
    print('Na: ', np.round(np.mean(na_c), 2))
if sum(zr_c) != 0:
    print('Zr: ', np.round(np.mean(zr_c), 2))
if sum(o_c) != 0:
    print('O: ', np.round(np.mean(o_c), 2))
if sum(p_c) != 0:
    print('P: ', np.round(np.mean(p_c), 2))
if sum(si_c) != 0:
    print('Si: ', np.round(np.mean(si_c), 2))
print('Total: ', np.round(total_c, 2))
print('------------------------------------------------------------------------------')

print('Writing data to the charge.dat file...')
   
filename = 'charge.dat'

with open(filename, 'w') as f:
    f.write('----------------------------Charge Analysis (Average)--------------------------------------\n')
    if sum(na_c) != 0:
        f.write('Na: ' + str(np.round(np.mean(na_c), 2)) + '\n')
    if sum(zr_c) != 0:
        f.write('Zr: ' + str(np.round(np.mean(zr_c), 2)) + '\n')
    if sum(o_c) != 0:
        f.write('O: ' + str(np.round(np.mean(o_c), 2)) + '\n')
    if sum(p_c) != 0:
        f.write('P: ' + str(np.round(np.mean(p_c), 2)) + '\n')
    if sum(si_c) != 0:
        f.write('Si: ' + str(np.round(np.mean(si_c), 2)) + '\n')
    f.write('Total: ' + str(np.round(total_c, 2)) + '\n')
    f.write('------------------------------------------------------------------------------------------\n')

print('------------------------------------------------------------------------------')

        
        







