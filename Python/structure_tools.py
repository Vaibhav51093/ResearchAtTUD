#######################################################
# python library for manipulating atoms
# requires: python-ase, python-asap, python-numpy
#
#######################################################

import math
import numpy as np
#import asap3 as asap
print("WARNING: asap3 not correctly installed for python3")
import ase
import ase.visualize
import ase.build
import ase.data 
import ase.neighborlist
import collections
import time
import os
from operator import itemgetter
#import spglib
from ase.parallel import paropen, parprint
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
import scipy.optimize
import spglib
import warnings
warnings.filterwarnings("ignore")





def symmetrise_cell(atoms, cell_tol=None):
    def_tol = {"lengths":0.001, "angles": 0.1}
    if cell_tol:
        for key in cell_tol:
            if key in def_tol:
                def_tol[key] = cell_tol[key]
            else:
                parprint("Warning: %s is not a valid key for cell_tol. Ignoring ..." %key)
    cell_tol = def_tol
    
    # before doing anything:
    # rotate a on x axis; the rotate around x axis such that b has not z component; now we are happy
    z = np.array([0, 0, 1])
    a, b, c = atoms.get_cell()
    atoms.rotate(a, "x", rotate_cell=True)
    theta = np.arccos( b.dot(z) / np.linalg.norm(z) / np.linalg.norm(z) )
    atoms.rotate(-theta, "x", rotate_cell=True)
    a, b, c = atoms.get_cell()
    if b.dot(z) < 0:
        atoms.rotate(180, "x", rotate_cell=True)        
    
    a0, b0, c0, alpha0, beta0, gamma0 = atoms.get_cell_lengths_and_angles()
    a, b, c = atoms.get_cell()
    if cell_tol["lengths"]:
        if 2*abs(a0 - b0)/(a0+b0) < cell_tol["lengths"]:
            if 2*abs(a0 - c0)/(a0+c0) < cell_tol["lengths"] or 2*abs(b0 - c0)/(b0+c0) < cell_tol["lengths"]:
                a_av = (a0 + b0 + c0)/3
                a = a_av / a0 * a
                b = a_av / b0 * b
                c = a_av / c0 * c
            else:
                a_av = (a0 + b0)/2
                a = a_av / a0 * a
                b = a_av / b0 * b

        elif 2*abs(a0 - c0)/(a0+c0) < cell_tol["lengths"]:
            a_av = (a0 + c0)/2
            a = a_av / a0 * a
            c = a_av / c0 * c

        elif 2*abs(b0 - c0)/(b0+c0) < cell_tol["lengths"]:
            a_av = (b0 + c0)/2
            b = a_av / b0 * b
            c = a_av / c0 * c
        atoms.set_cell([a, b, c], scale_atoms=True)

    a0, b0, c0, alpha0, beta0, gamma0 = atoms.get_cell_lengths_and_angles()
    if cell_tol["angles"]:
        if 2*abs(alpha0 - beta0) < cell_tol["angles"]:
            if 2*abs(alpha0 - gamma0) < cell_tol["angles"] or 2*abs(beta0 - gamma0) < cell_tol["angles"]:
                alpha0 = (alpha0 + beta0 + gamma0)/3
                beta0  = alpha0
                gamma0 = alpha0
            else:
                alpha0 = (alpha0 + beta0)/2
                beta0  = alpha0
        elif 2*abs(alpha0 - gamma0) < cell_tol["angles"]:
            alpha0 = (alpha0 + gamma0)/2
            gamma0 = alpha0
        elif 2*abs(beta0 - gamma0) < cell_tol["angles"]:
            beta0  = (beta0 + gamma0)/2
            gamma0 = beta0

    b1 = np.cos(gamma0/360 * 2*np.pi)*b0
    b2 = np.sqrt(b0**2 - b1**2)

    c1 = np.cos(beta0/360 * 2*np.pi)*c0
    c2 = (np.cos(alpha0/360 * 2*np.pi)*b0*c0 - b1*c1)/b2
    c3 = np.sqrt(c0**2 - c1**2 - c2**2)

    a = [a0, 0, 0]
    b = [b1, b2, 0]
    c = [c1, c2, c3]
    atoms.set_cell([a, b, c], scale_atoms=True)

    return atoms

    



                
                
            



def strain_tensor(eps=None):
    I   = np.diag([1, 1, 1])

    if type(eps) in [float, int]:
        D = np.diag([1, 1, 1])*eps/100.
    elif type(eps) in [np.ndarray, list]:
        D = np.zeros((3, 3))
        D[0][0] = eps[0]/100.
        D[1][1] = eps[1]/100.
        D[2][2] = eps[2]/100.
        D[1][2] = eps[3]/200.
        D[2][1] = eps[3]/200.
        D[0][2] = eps[4]/200.
        D[2][0] = eps[4]/200.
        D[0][1] = eps[5]/200.
        D[1][0] = eps[5]/200.
    else:
        print("Warning: eps must be a float, integer, list or np.array")
        print("The returned strain tensor is simply the unity matrix")
        return I

    return I+D


def deform(atoms, eps=1):
    ats  = atoms.copy()
    cell = atoms.get_cell()
    T    = strain_tensor(eps=eps)
    cell = T.dot(cell.T)
    ats.set_cell(cell.T, scale_atoms=True)
    return ats


def identify_deformation(atoms, atoms0):
    cell0 = atoms0.get_cell()
    cell0 = cell0.T
    cell = atoms.get_cell()
    cell = cell.T
    
    T = cell.dot(np.linalg.inv(cell0))
    D = T-np.diag([1, 1, 1])
    eps = [D[0][0], D[1][1], D[2][2], D[2][1]+D[1][2], D[2][0]+D[0][2], D[1][0]+D[0][1]]
    eps = 100*np.array(eps)
    eps = np.round(eps, 1)
    
    return eps



def count_species(atoms):
    Nx = {}
    for at in atoms:
        if at.symbol in Nx.keys():
            Nx[at.symbol] += 1
        else:
            Nx[at.symbol] = 1
    return Nx


def get_composition(atoms, normalize=False, string=False, fu=False):
    comp =  count_species(atoms)
    if fu == True:
        normalize = True
    if normalize == False:
        if string == True:
            keys = list(comp.keys())
            keys.sort()
            txt = ""
            for k in keys:
                txt += k
                txt += "%i_" %comp[k]
            txt = txt[:-1]
            return txt
        else:
            return comp
    else:
        species = comp.keys()
        natoms  = np.fromiter(comp.values(), dtype=int)
        if fu == True:
            nmin = GreatedCommonDivisor(natoms)
        else:
            nmin = natoms.min()
        natoms //= nmin
                
        comp = dict(zip(species, natoms))
        if string == True:
            keys = comp.keys()
            keys.sort()
            txt = ""
            for k in keys:
                txt += k
                txt += "%.3f_" %np.round(comp[k], 3)
                txt = txt[:-1]
            return txt
        else:
            return comp
        

def GreatedCommonDivisor(numbers):
    num = list(numbers[1:])
    
    gcd = numbers[0]
    while len(num) > 0:
        gcd = math.gcd(gcd, num[0])
        num = num[1:]
    return gcd
    
    
    
        

def remove(atoms, what=None):
    if type(what) == list:
        kill = what[:]
    if type(what) == str or type(what[0]) == str:
        if type(what) == str:
            what = [what]
        kill = []
        for at in atoms:
            if at.symbol in what:
                kill.append(at.index)
    else:
        print("what must be a list of integers, a string (chem element) or a list of strings ...")
        print("doing nothing")
        return atoms
                
    kill.sort()
    kill.reverse()
    for k in kill:
        del atoms[k]

    return atoms
                                                                        




def get_bond_lengths(atoms, cutoff=2.3):
    nl = get_ase_neighborlist(atoms, cutoff=cutoff)
    bonds = {}
    for at in atoms:
        i = at.index
        for j in nl[i]["i"]:
            k = nl[i]["i"].index(j)
            bond = [at.symbol, atoms[j].symbol]
            bond.sort()
            bond = bond[0]+"-"+bond[1]
            if bond not in bonds.keys():
                bonds[bond] = {"values":[] }
            bonds[bond]["values"].append(nl[i]["d"][k])
            
    # mean, stdev, maxdev
    
    for bond in bonds.keys():
        bonds[bond]["values"] = np.array(bonds[bond]["values"])
        bonds[bond]["mean"]   = bonds[bond]["values"].mean()
        bonds[bond]["stdev"]  = bonds[bond]["values"].std()
        bonds[bond]["maxdev"] = max([bonds[bond]["values"].mean()-bonds[bond]["values"].min(), 
                                     -bonds[bond]["values"].mean()+bonds[bond]["values"].max()])
    
    return bonds



def get_optimum_site(atoms, cutoff=3, p = (0, 0, 0)):
    
    def minimum_distance(p, atoms=None, cutoff=3):
        atoms.append(ase.Atom("H", p))
        nl = get_ase_neighborlist(atoms, cutoff=3)
        d  = min(nl[len(atoms)-1]["d"])
        del atoms[-1]
        return -d
    
    opt = scipy.optimize.minimize(minimum_distance, p, args= (atoms.copy(), cutoff))
    return opt
    


def get_coordination_numbers(atoms, cutoff = None):
    i = ase.neighborlist.neighbor_list('i', atoms, cutoff)
    coord = np.bincount(i)
    return coord

def get_rdf(atoms, cutoff = None, nbins=100):
    if not cutoff:
        a, b, c, alp, bet, gam = atoms.get_cell_lengths_and_angles()
        cutoff  = min([a, b, c])
        cutoff /= 2

    d = ase.neighborlist.neighbor_list('d', atoms, cutoff)
    h, bin_edges = np.histogram(d, bins=nbins)
    pdf = h/(4*np.pi/3*(bin_edges[1:]**3 - bin_edges[:-1]**3)) * atoms.get_volume()/len(atoms)
    r = (bin_edges[:-1] + bin_edges[1:])/2
    return r, pdf


def get_falsely_coordinated_atoms(atoms, cutoff=None, target = None):
    coord  = get_coordination_numbers(atoms, cutoff=cutoff)
    falsely = {}
    for key in target:
        falsely[key] = []

    for at in atoms:
        X = at.symbol
        i = at.index
        diff = coord[i] - target[X]
        if diff != 0:
            falsely[X].append([i, diff])

    return falsely


def get_local_bonding_environments(atoms, cutoff=None):
    NL = get_ase_neighborlist( atoms, cutoff=cutoff)
    bonds = []
    for i in range(len(atoms)):
        bonds.append({})
        X = atoms[i].symbol
        nl = NL[i]
        for j in nl["i"]:
            Y = atoms[j].symbol
            bond = [X, Y]
            bond.sort()
            bond = bond[0]+"-"+bond[1]
            if bond not in bonds[i].keys():
                bonds[i][bond] = 0
            bonds[i][bond] += 1

    return bonds
    


def count_falsely_coordinated_atoms(atoms, cutoff=None, target = None):
    falsely = get_falsely_coordinated_atoms(atoms, cutoff=cutoff, target = target)

    nfalsely = {"all":0}
    for key in falsely.keys():
        nfalsely[key] = len(falsely[key])
        nfalsely["all"] += nfalsely[key]
    
    return nfalsely








def get_asap_coordination_numbers(atoms, cutoff = None, bonds = None):
    if cutoff == None:
        print("Error, please specify a cutoff")
        return []

    if bonds != None:
        if type(bonds) != dict:
            print("Error, bonds need to be None or a dictionary")
            return []
        for key in bonds.keys():
            if type(bonds[key]) != list:
                bonds[key] = [bonds[key]]
            
        
    
    NL = asap.FullNeighborList(cutoff, atoms)
    NC = []
    if bonds == None:
        for nl in NL:
            NC.append(len(nl))
    else:
        for at in atoms:
            X = at.symbol
            i = at.index
            if X in bonds.keys():
                nc = 0
                for j in NL[i]:
                    if atoms[j].symbol in bonds[X]:
                        nc += 1
                NC.append(nc)
            else:
                NC.append(None)
    return NC
            


def get_undercoordinated_atoms(atoms, cn=None, cutoff=None):
    species = []
    for at in atoms:
        if at.symbol not in species:
            species.append(at.symbol)

    if type(cn) == int:
        tmp = {}
        for X in species:
            tmp[X] = cn
        cn = tmp.copy()
        

    if type(cutoff) in [int, float]:
        neighborlist = asap.FullNeighborList(cutoff, atoms)
        nblist = {}
        tmp    = {}
        for X in species:
            nblist[X] = neighborlist
            tmp[X] = cutoff
        cutoff = tmp.copy()
        
    elif type(cutoff) == dict:
        nblist = {}
        for X in species:
            nblist[X] = asap.FullNeighborList(cutoff[X], atoms)
        

    i_false = []
    for i in range(len(atoms)):
        X = atoms[i].symbol
        neighbors = nblist[X][i]
        if len(neighbors) < cn[X]:
            i_false.append(i)

    return i_false





def get_top_atom(atoms, axis="z"):
    atoms.wrap()
    z2i  = {"x":0, "y":1, "z":2}
    if type(axis) == str:
        axis = z2i[axis]
        z0 = 0
        i0 = 0
        pos = atoms.get_positions()
    else:
        z0 = -0.1
        i0 = 0
        pos = atoms.get_scaled_positions()
                        
    for i in range(len(pos)):
        if pos[i][axis] > z0:
            z0 = pos[i][axis]
            i0 = i
    return i0


def get_bottom_atom(atoms, axis="z"):
    atoms.wrap()
    z2i  = {"x":0, "y":1, "z":2}
    if type(axis) == str:
        axis = z2i[axis]
        z0 = np.linalg.norm(atoms.get_cell()[axis])
        i0 = 0
        pos = atoms.get_positions()
    else:
        z0 = 1.1
        i0 = 0
        pos = atoms.get_scaled_positions()
        
    for i in range(len(pos)):
        if pos[i][axis] < z0:
            z0 = pos[i][axis]
            i0 = i
    return i0

def get_bottom_atoms(atoms, axis="z", tol=0.5):
    z2i  = {"x":0, "y":1, "z":2}
    i0   = get_bottom_atom(atoms, axis=axis)
    if type(axis) == str:
        axis = z2i[axis]
        pos  = atoms.get_positions()
    else:
        cell = atoms.get_cell()
        pos  = atoms.get_scaled_positions()
        tol  = tol/np.linalg.norm(cell[axis])

    ibot = []
    for i in range(len(pos)):
        if abs(pos[i][axis] - pos[i0][axis]) < tol:
            ibot.append(i)
    return ibot

def get_top_atoms(atoms, axis="z", tol=0.5):
    z2i  = {"x":0, "y":1, "z":2}
    i0   = get_top_atom(atoms, axis=axis)
    if type(axis) == str:
        axis = z2i[axis]
        pos  = atoms.get_positions()
    else:
        cell = atoms.get_cell()
        pos  = atoms.get_scaled_positions()
        tol  = tol/np.linalg.norm(cell[axis])

    itop = []
    for i in range(len(pos)):
        if abs(pos[i][axis] - pos[i0][axis]) < tol:
            itop.append(i)
                                                                    
    return itop
                                            

def sort(atoms, order="312", reverse=False, pbc=None):
    return sort_atoms(atoms, order=order, reverse=reverse, pbc=pbc)


def sort_atoms(ats, order="zxy", reverse=False, pbc=None):
    atoms   = ats.copy()
    pbc_old = atoms.get_pbc()
    if pbc != None:
        atoms.set_pbc(pbc)

    tags = []
    if type(order) == int:
        for i in range(len(atoms)):
            d = atoms.get_distance(order, i, mic=True)
            if reverse == True:
                if d > 0.01:
                    d = 1/d
                else:
                    d = 9999999999
            tags.append(d)

    elif "m" in order:
        tags = atoms.get_masses()
        if reverse == True:
            tags = 1/tags

    elif "a" in order:
        tags = atoms.get_chemical_symbols()

    elif order == "tags":
        tags = atoms.get_tags()
        if reverse == True:
            for i in range(len(atoms)):
                if tags[i] > 0.001:
                    tags[i] = 1./tags[i]
                else:
                    tags[i] = 9999999999999999999999


    elif type(order) == str:
        tags = []
        pos  = atoms.get_positions()
        scp  = atoms.get_scaled_positions()
        for i in range(len(atoms)):
            tag = []
            for j in order:
                if j == "1":
                    t = scp[i][0]
                elif j == "2":
                    t = scp[i][1]
                elif j == "3":
                    t = scp[i][2]
                elif j == "x":
                    t = pos[i][0]
                elif j == "y":
                    t = pos[i][1]
                elif j == "z":
                    t = pos[i][2]
                if reverse == True:
                    t = 1/t
                elif type(reverse) == str:
                    if j in reverse:
                        if t > 0.01:
                            t = 1/t
                        else:
                            t = 999999999999
                               
                tag.append(t)
            tags.append(tag)
    elif type(order) == np.ndarray and len(order) == 3:
        atoms.append(ase.Atom("H", order))
        for i in range(len(atoms)-1):
            d =atoms.get_distance(-1, i, mic=True)
            if reverse == True:
                if d > 0.001:
                    d = 1./d
                else:
                    d = 9999999999999999999999
            tags.append(d)
        del atoms[-1]

        
    elif type(order) == list:
        tags = []        
        for at in atoms:
            tags.append(order.index(at.symbol))
        if reverse == True:
            print("Warning: the option reverse has no effect when chemical species are specified")
    
    atoms = ase.build.sort(atoms, tags=tags)
    atoms.set_pbc(pbc_old)
    return atoms
    


def reflect(ats, axis=None, center="origin"):
    if axis == None:
        return point_mirror(ats,  center=center)
    else:
        return plane_mirror(ats, axis=axis, center=center)



def point_mirror(ats,  center="origin"):
    atoms = ats.copy()
    a, b, c = atoms.get_cell()
    if type(center) == str:
        if center == "origin":
            center = np.array([0, 0, 0])
        if center == "coc":
            center = 0.5*(a+b+c)
    else:
        center = np.array(center)

    atoms.translate(-center)
    for at in atoms:
        at.position = -at.position
    atoms.translate(center)

    atoms.wrap(pbc=True)
    return atoms
            
    
def plane_mirror(ats, axis="z", center="origin"):
    atoms   = ats.copy()
    z2i     = {"x":0, "y":1, "z":2}
    axis    = z2i[axis]
    
    a, b, c = atoms.get_cell()
    if type(center) == str:
        if center == "origin":
            center = np.array([0, 0, 0])
        if center == "coc":
            center = 0.5*(a+b+c)
    else:
        center = np.array(center)

    
    atoms.translate(-center)
    for at in atoms:
        at.position[axis] = -at.position[axis]
    atoms.translate(center)

    atoms.wrap(pbc=True)
    return atoms



def get_ase_neighborlist(atoms, cutoff=3, distances=False):
    ind1, ind2, dist, vecs = ase.neighborlist.neighbor_list("ijdD", atoms, cutoff)
    NL = {}
    if ind1 != []:
        nl = [[ind1[0], ind2[0], dist[0], vecs[0]]]
        for i in range(1, len(ind1)):
            if ind1[i] == ind1[i-1]:
                nl.append([ind1[i], ind2[i], dist[i], vecs[i]])
            else:
                i1, i2, dis, vec = zip(*nl)
                NL[ind1[i-1]] = {"i": i2, "d":dis, "v": vec}
                nl = []
                nl.append([ind1[i], ind2[i], dist[i], vecs[i]])
                i1, i2, dis, vec = zip(*nl)

        i1, i2, dis, vec = zip(*nl)
        NL[ind1[i-1]] = {"i": i2, "d":dis, "v": vec}

    for i in range(len(atoms)):
        if i not in NL.keys():
            NL[i] = {"i": [], "d":[], "v": []}
            
    if distances == True:
        return NL, dist
    else:
        return NL




def get_ase_neighborshells(atoms, cutoff=None, NL=None, tol=0.2):
    NS = {} # NS[i] = {1:{"rav": rav, "d": [dk, dl, ... ], "i": [k,l, m, n], "atoms": ase.Atoms},
            #          2:{"rav": rav, "d": [dk, dl, ... ], "i": [k,l, m, n]}, ...
    if NL == None:
        NL = get_ase_neighborlist(atoms, cutoff=cutoff)

    for i in range(len(atoms)):
        shell = {}
        ind  = NL[i]["i"]
        dist = NL[i]["d"]
        vec  = NL[i]["v"]
        data = zip(ind, dist, vec)
        data.sort(key=itemgetter(1))
        
        nS = 1
        rS  = None
        for dat in data:
            j, d, v = dat
            at      = atoms[j]
            if rS == None:
                nn = 1   # number of atoms in shell
                rS = d   # av shell radius
                subshell = {"dav":0, "d":[d], "i":[j], "v":[v], "atoms":ase.Atoms([at])}
            
            elif abs(d-rS) > tol:
                subshell["dav"] = rS
                shell[nS] = subshell.copy()
                
                nS     += 1
                rS       = d
                nn       = 1
                subshell = {"dav":0, "d":[d], "i":[j], "v":[v], "atoms":ase.Atoms([at])}

            else:
                rS  = rS*nn + d
                nn  += 1
                rS /= nn
                subshell["i"].append(j)
                subshell["d"].append(d)
                subshell["v"].append(v)
                subshell["atoms"].append(at)

        subshell["dav"] = rS
        shell[nS] = subshell.copy()
                        
        
        NS[i] = shell.copy()
            
    return NS


def get_ase_rdf(atoms, cutoff=None, distances=None, bins=501, zeropad=True):
    if cutoff == None:
        a, b, c = atoms.get_cell()
        a0 = np.linalg.norm(a)
        b0 = np.linalg.norm(b)
        c0 = np.linalg.norm(c)
        cutoff = min([a0, b0, c0])
    if distances == None:
        dist   = ase.neighborlist.neighbor_list("d", atoms, cutoff)
        
    h, bin_edges = np.histogram(dist, bins=bins)
    x   = np.linspace(bin_edges[0], bin_edges[-1], bins)
    N   = len(atoms)
    rho = N/atoms.get_volume()
    pdf = 1./N/rho/(4/3.*np.pi)/(bin_edges[1:]**3 - bin_edges[:-1]**3) * h
    dx  = x[1]-x[0]

    if zeropad == True:
        N   = x[0]/dx
        x0  = np.linspace(0, x[0], N, endpoint=False)
        y0  = np.zeros(N)
        x   = np.array(list(x0)+list(x))
        pdf = np.array(list(y0)+list(pdf))
    
    
    return x, pdf



def get_asap_neighborlist(atoms, index=None,  cutoff=None, tol=0.33, nshell=2):
    if cutoff == None:
        a, b, c = atoms.get_cell()
        d = 0.5*(a+b+c)
        cutoff = np.linalg.norm(d)*0.99

    if index == None:
        index = range(len(atoms))
    elif type(index) == int:
        index = [index]

    # create neighborlist using asap
    ts = time.time()
    nblist = asap.FullNeighborList(cutoff, atoms=atoms)
    te = time.time()
    #print "ASAP time: %s" %(te-ts)
    
    # organize into new format:
    # neighbors[i][s]["radius"]  = shell radius
    # neighbors[i][s]["indices"] = list of neigbor indices of atom i in neighbor shell s
    # neighbors[i][s]["distances"] = list of corresponding distances of the neigbors  of atom i in shell s
    # neighbors[i][s]["vectors"] = same for vectors
    neighbors = {}
    for i in index:
        indices, vectors, distances = nblist.get_neighbors(i) #, cutoff)
        distances = list(np.sqrt(distances))
        indices   = list(indices)
        data = zip(indices, vectors, distances)
        data.sort(key=itemgetter(2))

        shells      = {}
        s_indices   = []
        s_distances = []
        s_vectors   = []
        rs     = data[0][-1]
        s      = 1
        for dat in data:
            index, vector, distance = dat
            if abs(distance-rs) > tol:
                shells[s]   = {"radius":rs, "indices":s_indices, "vectors": s_vectors, "distances": s_distances}
                s_indices   = [index]
                s_distances = [distance]
                s_vectors   = [vector]
                rs          = distance
                s          += 1
                if s > nshell:
                    break
            else:
                rs = distance
                s_indices.append(index)
                s_distances.append(distance)
                s_vectors.append(vector)

        neighbors[i] = shells
    return neighbors
                


def cut(ats, a=None, b=None, c=None, tol=0.01):
    ats.wrap(pbc=True)
    a0, b0 , c0 = ats.get_cell()
    V0 = ats.get_volume()
    N0 = len(ats)
    if a == None:
        a = a0.copy()
    if b == None:
        b = b0.copy()
    if c == None:
        c = c0.copy()

    if np.cross(a, b).dot(c) < 0:
        cell = [b, a, c]
    else:
        cell = [a, b, c]
    a, b, c = cell

    # basis transformation matrices from canonical cartesian basis
    BT = np.array([a, b, c]).T
    CT = np.linalg.inv(BT)
    
    car_positions = ats.get_positions()
    symbols       = ats.get_chemical_symbols()

    # scaled positions in new cell
    sca_positions = CT.dot(car_positions.T).T

    # fold all positions back into unit cell
    for i in range(len(sca_positions)):
        for j in range(3):
            if np.ceil(sca_positions[i][j]) > 1.01:
                sca_positions[i][j] -= np.floor(sca_positions[i][j])
            elif np.floor(sca_positions[i][j]) < - 0.01:
                sca_positions[i][j] -= np.floor(sca_positions[i][j])
            if sca_positions[i][j] > 0.9999:
                sca_positions[i][j] = 0
                                    
    # create new atoms object
    atoms = ase.Atoms([])
    for i in range(len(symbols)):    
        at = ase.Atom(symbols[i], sca_positions[i])
        atoms.append(at)
    atoms.set_cell([a, b, c], scale_atoms=True)
    atoms.set_pbc(True)
    
    # remove duplicates
    kill = []
    for i in range(len(atoms)-1):
        if i not in kill:
            for j in range(i+1, len(atoms)):
                if j not in kill:
                    if atoms.get_distance(i, j, mic=True) < 0.1:
                        kill.append(j)
    kill.sort()
    kill.reverse()
    for k in kill:
        del atoms[k]
        
    # check whether N/V = N0/V0
    V = atoms.get_volume()
    N = len(atoms)
    while abs(N/V - N0/V0) > 0.001:
        # add old basis vectors +/- to existing positions
        # and keep if resulting position is still inside the cell!
        for n in range(N):
            at = atoms[n]
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        if abs(i)+abs(j)+abs(k) != 0:
                            atoms.append(at)
                            atoms[-1].position += i*a0 + j*b0 + k*c0
                            pos  = atoms[-1].position.copy()
                            scp  = CT.dot(pos)
                            for l in range(3):
                                if scp[l] > 0.9999 or scp[-l] < -0.0001:
                                    del atoms[-1]
                                    break
                            
        # remove duplicates again
        kill = []
        for i in range(len(atoms)-1):
            if i not in kill:
                for j in range(i+1, len(atoms)):
                    if j not in kill:
                        if atoms.get_distance(i, j, mic=True) < 0.1:
                            kill.append(j)
        kill.sort()
        kill.reverse()
        for k in kill:
            del atoms[k]

        V = atoms.get_volume()
        N = len(atoms)

    atoms.set_pbc(True)
    return atoms



    
def get_csl(ats, nxn=[17, 17], axis="z", max_shells=12):
    # keep only non-equivalent atoms in the cell

    # keep only the outermost layer
    ibot    = get_bottom_atoms(ats, axis=axis) 
    a, b, c = ats.get_cell()
    atoms   = ase.Atoms([])
    for i in ibot:
        atoms.append(ats[i])
    atoms.set_cell([a, b, c])
    atoms.translate(-np.array([0, 0, atoms[0].position[2]]))
    a, b, c = atoms.get_cell()
    A0 = np.linalg.norm(np.cross(a, b))

    # create supercell by in-plane repetition
    atoms = atoms.repeat(nxn+[1])

    # get atoms closest to the plane center
    a, b, c = atoms.get_cell()
    atoms.set_pbc(False)
    com = atoms.get_center_of_mass()
    atoms.append(ase.Atom("H", com))
    d0 = atoms.get_distance(0, -1)
    icenter = 0
    for i in range(1, len(atoms)-1):
        d = atoms.get_distance(i, -1)
        if d < d0:
            d0 = d
            icenter = i
   
    # move this atom to origin (center atoms object at zero)
    atoms.translate(-atoms[icenter].position.copy())
    del atoms[-1]
    
    # compute ditances between all atoms and the atom at the origin
    tmp = []
    for at in atoms:
        d = np.linalg.norm(at.position)
        i = at.index
        tmp.append([i, d])
    tmp.sort(key=itemgetter(1))
    del tmp[0]

    # organize into shells
    nshells = []
    shell = []
    d0 = tmp[0][1]
    elem = 1
    for t in tmp:
        i = t[0]
        d = t[1]
        if abs(d-d0) < 0.1:
            shell.append(i)
            atoms[i].symbol = ase.data.chemical_symbols[elem]
        else:
            elem += 1
            atoms[i].symbol = ase.data.chemical_symbols[elem]
            nshells.append(shell)
            shell = [i]
            d0    = t[1]
    nshells.append(shell)
    ase.visualize.view(atoms)

    # compute rotation angles that lead to a coincidence
    # also compute the corresponding sigma value
    # 1. construct a right-handed cell by combining two vectors linear independent vectors (v1 from shell[i], v2 from shell[j])
    # 2. compute all angles between v1 and all other vectors within the same shell[i], theta_1
    # 3. compute all angles between v2 and all other vectors within the same shell[j], theta_2
    # 4. if there are common angles, theta_1 == theta_2 == theta, this rotation yields a coincidence!
    # 5. the sigma value is A/A0

    # symmetry operations on first neighbour shell
    c = np.cross(a, b)
    c = c/np.linalg.norm(c)    
    theta_sym = [0]
    i  = nshells[0][0]
    v1 = atoms[i].position.copy()
    theta1 = []
    for j in nshells[0]:
        # rotated cell vector 1, u1
        u1  = atoms[j].position.copy()
        cos_theta = round(v1.dot(u1)/v1.dot(v1), 5)
        theta = np.round(360 * np.arccos(cos_theta)/2/np.pi, 5)
        if abs(theta-90) < 0.01:
            print(c)
            if np.cross(v1, u1).dot(c) > 0:
                theta = 90.
            else:
                theta = 270.
        theta_sym.append(theta)
    #print theta_sym

    # find csl rotation angles ...
    csl = {}
    for shell1 in nshells[:max_shells]:
        # cell vector 1, v1
        i  = shell1[0]
        v1 = atoms[i].position.copy()
        theta1  = []
        v1_rot  = []
        for j in shell1:
            # rotated cell vector 1, u1
            u1  = atoms[j].position.copy()
            cos_theta = round(v1.dot(u1)/v1.dot(v1), 5)
            theta     = np.round(360 * np.arccos(cos_theta)/2/np.pi, 5)
            if abs(theta-90) < 0.01:
                if np.cross(v1, u1).dot(c) > 0:
                    theta =90.
                else:
                    theta =270.
                                
            if theta not in theta_sym:
                theta1.append(theta)
                v1_rot.append(u1)

        # cell vector 2, v2
        for shell2 in nshells[:max_shells]:
            for i in shell2:
                v2        = atoms[i].position.copy()
                # check that angle of cell is between 60 and 120 degrees
                cos_alpha = round(v1.dot(v2)/np.linalg.norm(v1)/np.linalg.norm(v2), 5)
                alpha     = 360/(2*np.pi)*np.arccos(cos_alpha)
                if abs(alpha-90) < 0.01:
                    if np.cross(v1, u1).dot(c) > 0:
                        alpha = 90.
                    else:
                        alpha = 270.
                                    
                if 59.999 < alpha < 120.001:
                    for j in shell2:
                        # rotated cell vector 2, u2
                        u2 = atoms[j].position.copy()
                        cos_theta = round(v2.dot(u2)/v2.dot(v2), 5)
                        theta = np.round(360 * np.arccos(cos_theta)/2/np.pi, 5)
                        if abs(theta-90) < 0.01:
                            if np.cross(v1, u1).dot(c) > 0:
                                theta =90.
                            else:
                                theta =270.
                                                                
                        if theta in theta1 and ( (theta-alpha) > 0.1 or abs(np.linalg.norm(v1)-np.linalg.norm(v2)) > 0.1) :
                            k  = theta1.index(theta)
                            u1 = v1_rot[k]
                            cos_alpha_rot = round(u1.dot(u2)/np.linalg.norm(u1)/np.linalg.norm(u2), 5)
                            alpha_rot     = 360/(2*np.pi)*np.arccos(cos_alpha_rot)
                            if abs(alpha-alpha_rot) < 0.1:
                                A = round(np.linalg.norm(np.cross(v1, v2)), 0)
                                sigma = int(round(A/A0, 0))
                                if sigma % 2 == 1:
                                    if sigma not in csl.keys():
                                        csl[sigma] = [theta, alpha, v1, v2]                                    
                                    elif abs(alpha-90) < 0.1:
                                        csl[sigma] = [theta, alpha, v1, v2]
                                    elif abs(alpha-120) < 0.1 and abs(csl[sigma][1]-90) > 0.1:
                                        csl[sigma] = [theta, alpha, v1, v2]
                                    elif abs(alpha-60) < 0.1 and abs(csl[sigma][1]-90) > 0.1 and abs(csl[sigma][1]-120) > 0.1:
                                        csl[sigma] = [theta, alpha, v1, v2]
                        
                            

    csl = collections.OrderedDict(sorted(csl.items()))
    for sigma in csl.keys():
        theta, alpha, v1, v2 = csl[sigma]
        a1 = atoms.copy()
        a2 = atoms.copy()
            
        for at in a1:
            at.symbol ="Si"
            at.position[2]+=1
            
        for at in a2:
            at.symbol ="C"
            at.position[2]+=3
        a2.rotate(theta, "z")
        a = a1+a2
        a.set_cell([v1, v2, c])
        csl[sigma].append(a)
        

    return csl
    
    


def get_polygons(atoms, index=None, maxdepth=10, cutoff= 1.8, npoly=3):
    nblist = asap.FullNeighborList(cutoff, atoms)
    polygons = []

    if index == None:
        index = range(len(atoms))
    if type(index) == int:
        index = [index]
    
    
    for i0 in index:
        paths = [[i0]]
        # construct all closed path of max length 'maxdepth'
        closed_paths = []
        # two steps through neighbor lists, avoiding to get back after the second step!
        for d in range(2):
            tmp = []
            for path in paths:
                p0 = path[-1]
                for n in nblist[p0]:
                    if n not in path:
                        tmp.append(path+[n])
            paths = tmp[:]                                

        # continue until maxdepth reached
        for d in range(2, maxdepth):
            tmp = []
            for path in paths:
                p0 = path[-1]
                for n in nblist[p0]:
                    if n not in path[1:]:
                        if n == i0:
                            closed_paths.append(path+[n])
                        else:
                            tmp.append(path+[n])
                        
            paths = tmp[:]
            
        l0 = maxdepth+1
        p  = []
        for path in closed_paths:
            p.append(len(path)-1)
        p.sort()
        P = []
        for i in range(min(npoly, len(p)/2)):
            P.append(p[2*i])
        
        polygons.append([ i0, P])


    # return only non-equivalent paths
    p_sorted = []
    for p in closed_paths:
        p = p[:-1]
        p.sort()
        p = np.array(p)
        accept = True
        for q in p_sorted:
            if len(q) == len(p):
                if np.linalg.norm(q-p) < 0.01:
                    accept=False
        if accept == True:
            p_sorted.append(p)
                                                                            
    return polygons,  closed_paths, p_sorted
    
    
def get_vasp_dos(a=None, loc=None):
    if loc == None:
        loc = os.getcwd()
    try:
        f = open("%s/DOSCAR" %loc, "r")
    except:
        print("DOSCAR available in %s" %loc)

    for i in range(5):
        f.readline()
        
    header = f.readline()
    header = header.split()
    npts = int(header[2])
    fermilevel = float(header[3])

    dos = []
    if a == None:
        for i in range(npts):
            line = f.readline()
            line = line.strip().split()
            dos.append(map(float, line))

        dos = np.array(dos)
        e, dos, int_dos = dos.T
        e -= fermilevel
        return {"e":e, "t":dos, "int":int_dos}
    
    else:
        for i in range(npts+1):
            line = f.readline()
        
        for i in range(a*(npts+1)):
            line = f.readline()

        for i in range(npts):
            line = f.readline()
            line = line.strip().split()
            dos.append(map(float, line))
        
        dos = np.array(dos)
        e, s, px, py, pz, d1, d2, d3, d4, d5 = dos.T
        p = px + py +pz
        d = d1 + d2 + d3 + d4 + d5
        t = s + p + d
        e -= fermilevel
    
        return {"e":e, "t":t, "s":s, "p":p, "d":d}
                                                                                                                                    


def get_layers(atoms, layer_dir=3, tol=0.75):
    species = []
    for at in atoms:
        if at.symbol not in species:
            species.append(at.symbol)
            
    
    layers_pos = {}
    layers_indices = {}
    for X in species:
        layers_pos[X] = []
        layers_indices[X] = []
        
    for at in atoms:
        p = at.position[layer_dir-1]
        X = at.symbol
        
        if layers_pos[X] == []:
            layers_pos[X].append(p)
            layers_indices[X].append([at.index])

        else:
            added = False
            for q in layers_pos[X]:
                if abs(p-q) < tol:
                    i = layers_pos[X].index(q)
                    layers_indices[X][i].append(at.index)
                    qq = 0
                    for j in layers_indices[X][i]:
                        qq += atoms[j].position[layer_dir-1]
                    qq /= len(layers_indices[X][i])
                    layers_pos[X][i] = qq
                    added = True
                    break
            if added == False:
                layers_pos[X].append(p)
                layers_indices[X].append([at.index])

    positions = []
    indices   = []
    symbols   = []
    for X in species:
        positions += layers_pos[X]
        indices   += layers_indices[X]
        for i in range(len(layers_pos[X])):
            symbols.append(X)
        
    layers = list(zip(positions, symbols, indices))
    layers.sort(key=itemgetter(0))
    positions, symbols, indices = zip(*layers)

    return positions, symbols, indices
            
    


def intercalate(atoms, X="F", Y=["O"], nmax=8, stol=1e-2, atol=5e-1):
    comp = SymmetryEquivalenceCheck(stol=stol)
    if type(Y) != list:
        Y = [Y]
        
    composition = get_composition(atoms)
    if X not in composition.keys():
        composition[X] = 0
        
    if composition[X] >= nmax:
        return

    dataset = spglib.get_symmetry_dataset(atoms, symprec=stol, angle_tolerance=atol, hall_number=0)
    equivalent_atoms = dataset["equivalent_atoms"]

    checked = []
    for ea in equivalent_atoms:
        if atoms[ea].symbol in Y and ea not in checked:
            checked.append(ea)
            tmp = atoms.copy()
            tmp[ea].symbol = X
        
            content = os.listdir(".")
            accept  = True
            counter = 0
            for c in content:
                if "POSCAR.intercalated" in c:
                    ref = ase.io.read(c)
                    if comp.compare(tmp, ref) == True:
                        accept = False
                        break
                    else:
                        counter += 1
            if accept == True:
                tmp.translate(-tmp[-1].position)
                tmp.wrap(pbc=True)
                print("Writing structure no. %i" %counter)
                tmp.write("POSCAR.intercalated.%i.vasp" %counter, sort=True, vasp5=True)
                intercalate(tmp, X=X, Y=Y, nmax=nmax, stol=stol,  atol=atol)

    return



def deintercalate(atoms, X="Na", nmax=0, stol=1e-2, atol=1e-1):
    comp = SymmetryEquivalenceCheck(stol=stol)
    composition = get_composition(atoms)
    if X not in composition.keys():
        composition[X] = 0
    if type(nmax) == float:
        nmax = int(composition[X]*nmax)
    #print nmax, composition[X]
    if composition[X] <= nmax:
        return

    dataset = spglib.get_symmetry_dataset(atoms, symprec=stol, angle_tolerance=atol, hall_number=0)
    equivalent_atoms = dataset["equivalent_atoms"]

    checked = []
    for ea in equivalent_atoms:
        if atoms[ea].symbol == X and ea not in checked:
            checked.append(ea)
            tmp = atoms.copy()
            del tmp[ea]

            content = os.listdir(".")
            accept  = True
            counter = 0
            for c in content:
                if "POSCAR.deintercalated" in c:
                    ref = ase.io.read(c)
                    if comp.compare(tmp, ref) == True:
                        accept = False
                        break
                    else:
                        counter += 1
        
            if accept == True:
                tmp.translate(-tmp[-1].position)
                tmp.wrap(pbc=True)
                print("Writing structure no. %i" %counter)
                tmp.write("POSCAR.deintercalated.%i.vasp" %counter, sort=True, vasp5=True)
                deintercalate(tmp, X=X, nmax=nmax, stol=stol,  atol=atol)


    return
                                                                                                                                                                                                                                                                                                                        


def remove_vacancy_placeholder(atoms, X="X", rescale=True, maxdist=2.5, mindist=1.8):
    while X in atoms.get_chemical_symbols():
        species = atoms.get_chemical_symbols()
        i = species.index(X)
        del atoms[i]

            
    # rescale
    atoms.wrap(pbc=True)
    atoms.center()
    a, b, c = atoms.get_cell()
    atoms.rotate(a, "x", rotate_cell=True)
    a, b, c = atoms.get_cell()
    z   = np.array([0, 0, 1])
    phi = np.arccos(b.dot(z)/np.linalg.norm(b))* 360 / 2 /np.pi
    atoms.rotate(phi-90, "x", rotate_cell=True)

    a, b, c = atoms.get_cell()
    if c.dot(z) < 0:
        atoms.rotate(180, "x", rotate_cell=True)

    
    
    if rescale == True:
        for order in [3, 2, 1]:
            atoms = sort(atoms, order=str(order))
            positions, symbols, indices = get_layers(atoms, layer_dir=order)
            for i in range(len(positions)-1):
                layer_sep = positions[i+1] - positions[i]
                if layer_sep > maxdist:
                    cell = atoms.get_cell()
                    c0   = np.linalg.norm(cell[order-1])
                                      
                    i0 = min(indices[i+1])
                    for j in range(i0, len(atoms)):
                        atoms[j].position -= (layer_sep - mindist)*cell[order-1]/c0
            
                    
                    cell[order-1] = (c0 - (layer_sep - mindist)) * cell[order-1]/c0
                    atoms.set_cell(cell, scale_atoms = False)
            atoms.wrap(pbc=True)
                        
    return atoms
