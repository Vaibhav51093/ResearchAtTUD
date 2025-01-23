# K point convergence and parameters important for VASP bulk calculations 
from pyiron import Project, ase_to_pyiron
import matplotlib.pyplot as plt
import numpy as np
from pyiron import Project
from ase.io import read, write
from pyiron import ase_to_pyiron, pyiron_to_ase
import ase
import os
import time
import tarfile
from vasp_out import read_IBZKBT


class k_sampling():

    def k_points(struct, rang=range(1,21,1)):
        
        # Input is just structure from pyiron, how much division you want  

        k_points_test = []

        struct_1 = pyiron_to_ase(struct)
        cell_info = struct_1.get_cell_lengths_and_angles()
        lattice_para = cell_info[0:3]
        
        # Now the task is to descritize reciprocal lattice by dividing reci lattice 
        division = rang
        for k in division:
            if lattice_para[0]==lattice_para[1] and lattice_para[0]==lattice_para[2]:
                reci_space = (2*(22/7))/lattice_para[0]
                space = reci_space/k
                k_point = reci_space/space
                k_points_test.append('{} {} {}'.format(int(round(k_point)),int(round(k_point)),int(round(k_point))))
                #print('K-points = {} x {} x {}'.format(int(round(k_point)),int(round(k_point)),int(round(k_point))))        

            elif lattice_para[0]==lattice_para[1] and lattice_para[0] != lattice_para[2]:
                reci_space = (2*(22/7))/lattice_para[0]
                space = reci_space/k
                k_point = reci_space/space
                reci_space_2 = (2*(22/7))/lattice_para[2]
                space_2 = reci_space_2/k
                k_point_2 = reci_space_2/space
                k_points_test.append('{} {} {}'.format(int(round(k_point)),int(round(k_point)),int(round(k_point_2))))
                #print('K-points = {} x {} x {}'.format(int(round(k_point)),int(round(k_point)),int(round(k_point_2))))

            else:
                reci_space = (2*(22/7))/lattice_para[0]
                space = reci_space/k
                k_point = reci_space/space
                reci_space_1 = (2*(22/7))/lattice_para[1]
                space_1 = reci_space_1/k
                k_point_1 = reci_space_1/space
                reci_space_2 = (2*(22/7))/lattice_para[2]
                space_2 = reci_space_2/k
                k_point_2 = reci_space_2/space
                k_points_test.append('{} {} {}'.format(int(round(k_point)),int(round(k_point_1)),int(round(k_point_2))))


        return k_points_test

    # K_density per angstrom 
    def k_density(N_x, N_y, N_z, a, b, c):

        # Provide K points in x, y and z direction 
        
        den_x = ((N_x*a*7)/(2*22)) 
        den_y = ((N_y*b*7)/(2*22))
        den_z = ((N_z*c*7)/(2*22))
        return den_x,den_y,den_z

    # K_density to K_points 
    def k_points_density(dens, struct):
        
        # Provide k-point density and structure in pyiron 
        
        k_points_test = []
        struct_1 = pyiron_to_ase(struct)
        cell_info = struct_1.get_cell_lengths_and_angles()
        lattice_para = cell_info[0:3]
        a = lattice_para[0]
        b = lattice_para[1]
        c = lattice_para[2]
        N_x = ((dens*2*22)/(7*a))
        N_y = ((dens*2*22)/(7*b))
        N_z = ((dens*2*22)/(7*c))
        #print(N_x)
        k_points_test.append('{} {} {}'.format(int(round(N_x)),int(round(N_y)),int(round(N_z))))
        return k_points_test


class convergence():

    # While doing convergence test one can get values or just energy per atom 
    # 1. Tot energy/atom (eV)
    # 2. Force
    # 3. IBZKPT 
    # 4. energy difference if multiple inputs, set diff_1=true
    def enrg_force_ibz_diff(pyiron_job):

        # Get energy, force, irr. BZ points for each job 
        job = pyiron_job
        tot_energy_even = job['output/generic/energy_pot'][-1]/len(job.output.positions[-1])
        resultant_force = (np.sum(job.output.forces[-1]))  # force sum in x, y, z direction 
        max_force = job.output.force_max[-1]               # Max force in last iteration 
        ibz_2 = read_IBZKBT.IBZ_vasp(job)
            
        return tot_energy_even, resultant_force, max_force, ibz_2

    def convergence(pyiron_job,rang,diff_1=False):

        # for K-point and Ecut convergence seperately use this one
        # Input needed :- 1. pyiron job with placeholder 
        #                 2. Range of jobs 
        #                 3. want to calculate difference
 
        tot_energy_even = []
        force = []
        ibz = []
        f_max = []
        diff = [] 
        job = pyiron_job
        for i in rang:
            tot_energy_even_2,resultant_force_2,max_force_2,ibz_3 = enrg_force_ibz_diff(pyiron_job,ibz_1=False)
            tot_energy_even.append(tot_energy_even_2)
            force.append(resultant_force_2)
            ibz.append(ibz_3)
            f_max.append(max_force_2)
        if diff_1 == True:
            for i in range(len(tot_energy_even)):
                diff.append(abs(tot_energy_even[i]-tot_energy_even[len(tot_energy_even)-1]))
        return tot_energy_even, diff, ibz, force, f_max 

    def kpoint_ecut_convergence(pyiron_job,rang,range_2,diff_1=False):

        # for Ecut and  kpoint convergence together 
        # Input needed :- 1. pyiron job with two palceholder as per job name
        #                    rang = ecut, range_2 = kpoints or vice versa   
        #                 2. Range of jobs 
        #                 3. want to calculate difference

        tot_energy_even_n = []
        force_n = []
        ibz_n = []
        f_max_n = []
        diff_n = [] 
        job = pyiron_job

        for j in rang:

            tot_energy_even = []
            force = []
            ibz = []
            f_max = []
            diff = [] 

            for i in rang:
                tot_energy_even_2,resultant_force_2,max_force_2,ibz_3 = enrg_force_ibz_diff(pyiron_job,ibz_1=False)
                tot_energy_even.append(tot_energy_even_2)
                force.append(resultant_force_2)
                ibz.append(ibz_3)
                f_max.append(max_force_2)
            if diff_1 == True:
                for i in range(len(tot_energy_even)):
                    diff.append(abs(tot_energy_even[i]-tot_energy_even[len(tot_energy_even)-1]))

            return tot_energy_even, diff, ibz, force, f_max

        return tot_energy_even_n, diff_n, ibz_n, force_n, f_max_n

class vasp_parallel():
    # NCORE and KPAR 
    def ncore_kpar(core_to_us=84,ir_point=84):
        if (ir_point % 2) == 0:
            k_list = [1,2,3,4,5,6,7,8,9,10,11,12]
            n_list = [2,4,3,6,8,12]
            kpar_n = []
            ncore_n = []
            for i in k_list:
                kp = (ir_point)/(core_to_us/i)
                if (kp).is_integer() == True:
                    kpar=i
                    kpar_n.append(i)
                    ir_per_div = ir_point/kpar
                    for k in n_list:
                        if (ir_per_div/k).is_integer() == True:
                            n_core = k
                            ncore_n.append(k)
            return kpar_n, ncore_n 




# Modify and add other functions realted to     