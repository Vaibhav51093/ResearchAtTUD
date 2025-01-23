from scipy.optimize import curve_fit
#from pyiron import Project
#from pyiron import ase_to_pyiron, pyiron_to_ase
import shutil
import glob
from scipy import stats
from scipy.optimize import curve_fit
import numpy as np 
import scipy.constants as const
import pandas as pd
from scipy.optimize import curve_fit
import lammps_logfile
import numpy as np 
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *


def compute_msd_sigma(filee_list=['msd_na_nasi_523_1_0_5000.txt', 'msd_na_nasi_523_1_0_5000.txt', 'msd_na_nasi_523_1_0_5000.txt', 'msd_na_nasi_523_1_0_5000.txt', 'msd_na_nasi_523_1_0_5000.txt'], output_file='msd_glass_523k_avg_sigma.png', xlimit=[0, 2], ylimt=True, y_limit=[0, 900]):
    
    import numpy as np

    MSD_n_x_1 = []
    MSD_n_x_2 = []
    MSD_n_x_3 = []
    MSD_n_x_4 = []
    MSD_n_x_5 = []
    MSD_n_y_1 = []
    MSD_n_y_2 = []
    MSD_n_y_3 = []
    MSD_n_y_4 = []
    MSD_n_y_5 = []
    MSD_n_z_1 = []
    MSD_n_z_2 = []
    MSD_n_z_3 = []
    MSD_n_z_4 = []
    MSD_n_z_5 = []
    MSD_n_xy_1 = []
    MSD_n_xy_2 = []
    MSD_n_xy_3 = []
    MSD_n_xy_4 = []
    MSD_n_xy_5 = []
    MSD_n_xz_1 = []
    MSD_n_xz_2 = []
    MSD_n_xz_3 = []
    MSD_n_xz_4 = []
    MSD_n_xz_5 = []
    MSD_n_yz_1 = []
    MSD_n_yz_2 = []
    MSD_n_yz_3 = []
    MSD_n_yz_4 = []
    MSD_n_yz_5 = []
    time = []

    for i in range(1):
        if i == 0:

            file_data_1 = np.loadtxt(filee_list[0])
            file_data_2 = np.loadtxt(filee_list[1])
            file_data_3 = np.loadtxt(filee_list[2])
            file_data_4 = np.loadtxt(filee_list[3])
            file_data_5 = np.loadtxt(filee_list[4])
            
            t = file_data_1[:,0]
            time.append(t)
        
            netMSD_x_1 = file_data_1[:,7]
            netMSD_x_2 = file_data_2[:,7]
            netMSD_x_3 = file_data_3[:,7]
            netMSD_x_4 = file_data_4[:,7]
            netMSD_x_5 = file_data_5[:,7]

            netMSD_y_1 = file_data_1[:,8]
            netMSD_y_2 = file_data_2[:,8]
            netMSD_y_3 = file_data_3[:,8]
            netMSD_y_4 = file_data_4[:,8]
            netMSD_y_5 = file_data_5[:,8]

            netMSD_z_1 = file_data_1[:,9]
            netMSD_z_2 = file_data_2[:,9]
            netMSD_z_3 = file_data_3[:,9]
            netMSD_z_4 = file_data_4[:,9]
            netMSD_z_5 = file_data_5[:,9]

            netMSD_xy_1 = file_data_1[:,10]
            netMSD_xy_2 = file_data_2[:,10]
            netMSD_xy_3 = file_data_3[:,10]
            netMSD_xy_4 = file_data_4[:,10]
            netMSD_xy_5 = file_data_5[:,10]

            netMSD_xz_1 = file_data_1[:,11]
            netMSD_xz_2 = file_data_2[:,11]
            netMSD_xz_3 = file_data_3[:,11]
            netMSD_xz_4 = file_data_4[:,11]
            netMSD_xz_5 = file_data_5[:,11]

            netMSD_yz_1 = file_data_1[:,12]
            netMSD_yz_2 = file_data_2[:,12]
            netMSD_yz_3 = file_data_3[:,12]
            netMSD_yz_4 = file_data_4[:,12]
            netMSD_yz_5 = file_data_5[:,12]


            MSD_n_x_1.append(netMSD_x_1)
            MSD_n_x_2.append(netMSD_x_2)
            MSD_n_x_3.append(netMSD_x_3)
            MSD_n_x_4.append(netMSD_x_4)
            MSD_n_x_5.append(netMSD_x_5)

            MSD_n_y_1.append(netMSD_y_1)
            MSD_n_y_2.append(netMSD_y_2)
            MSD_n_y_3.append(netMSD_y_3)
            MSD_n_y_4.append(netMSD_y_4)
            MSD_n_y_5.append(netMSD_y_5)

            MSD_n_z_1.append(netMSD_z_1)
            MSD_n_z_2.append(netMSD_z_2)
            MSD_n_z_3.append(netMSD_z_3)
            MSD_n_z_4.append(netMSD_z_4)
            MSD_n_z_5.append(netMSD_z_5)

            MSD_n_xy_1.append(netMSD_xy_1)
            MSD_n_xy_2.append(netMSD_xy_2)
            MSD_n_xy_3.append(netMSD_xy_3)
            MSD_n_xy_4.append(netMSD_xy_4)
            MSD_n_xy_5.append(netMSD_xy_5)

            MSD_n_xz_1.append(netMSD_xz_1)
            MSD_n_xz_2.append(netMSD_xz_2)
            MSD_n_xz_3.append(netMSD_xz_3)
            MSD_n_xz_4.append(netMSD_xz_4)
            MSD_n_xz_5.append(netMSD_xz_5)

            MSD_n_yz_1.append(netMSD_yz_1)
            MSD_n_yz_2.append(netMSD_yz_2)
            MSD_n_yz_3.append(netMSD_yz_3)
            MSD_n_yz_4.append(netMSD_yz_4)
            MSD_n_yz_5.append(netMSD_yz_5)

    MSD_n_x_avg_1 = np.mean(MSD_n_x_1,axis=0)
    MSD_n_x_avg_2 = np.mean(MSD_n_x_2,axis=0)
    MSD_n_x_avg_3 = np.mean(MSD_n_x_3,axis=0)
    MSD_n_x_avg_4 = np.mean(MSD_n_x_4,axis=0)
    MSD_n_x_avg_5 = np.mean(MSD_n_x_5,axis=0)

    MSD_n_y_avg_1 = np.mean(MSD_n_y_1,axis=0)
    MSD_n_y_avg_2 = np.mean(MSD_n_y_2,axis=0)
    MSD_n_y_avg_3 = np.mean(MSD_n_y_3,axis=0)
    MSD_n_y_avg_4 = np.mean(MSD_n_y_4,axis=0)
    MSD_n_y_avg_5 = np.mean(MSD_n_y_5,axis=0)

    MSD_n_z_avg_1 = np.mean(MSD_n_z_1,axis=0)
    MSD_n_z_avg_2 = np.mean(MSD_n_z_2,axis=0)
    MSD_n_z_avg_3 = np.mean(MSD_n_z_3,axis=0)
    MSD_n_z_avg_4 = np.mean(MSD_n_z_4,axis=0)
    MSD_n_z_avg_5 = np.mean(MSD_n_z_5,axis=0)

    MSD_n_xy_avg_1 = np.mean(MSD_n_xy_1,axis=0)
    MSD_n_xy_avg_2 = np.mean(MSD_n_xy_2,axis=0)
    MSD_n_xy_avg_3 = np.mean(MSD_n_xy_3,axis=0)
    MSD_n_xy_avg_4 = np.mean(MSD_n_xy_4,axis=0)
    MSD_n_xy_avg_5 = np.mean(MSD_n_xy_5,axis=0)

    MSD_n_xz_avg_1 = np.mean(MSD_n_xz_1,axis=0)
    MSD_n_xz_avg_2 = np.mean(MSD_n_xz_2,axis=0)
    MSD_n_xz_avg_3 = np.mean(MSD_n_xz_3,axis=0)
    MSD_n_xz_avg_4 = np.mean(MSD_n_xz_4,axis=0)
    MSD_n_xz_avg_5 = np.mean(MSD_n_xz_5,axis=0)

    MSD_n_yz_avg_1 = np.mean(MSD_n_yz_1,axis=0)
    MSD_n_yz_avg_2 = np.mean(MSD_n_yz_2,axis=0)
    MSD_n_yz_avg_3 = np.mean(MSD_n_yz_3,axis=0)
    MSD_n_yz_avg_4 = np.mean(MSD_n_yz_4,axis=0)
    MSD_n_yz_avg_5 = np.mean(MSD_n_yz_5,axis=0)

    msd_n_x_avg = np.mean([MSD_n_x_avg_1,MSD_n_x_avg_2,MSD_n_x_avg_3,MSD_n_x_avg_4,MSD_n_x_avg_5], axis=0)
    msd_n_y_avg = np.mean([MSD_n_y_avg_1,MSD_n_y_avg_2,MSD_n_y_avg_3,MSD_n_y_avg_4,MSD_n_y_avg_5], axis=0)
    msd_n_z_avg = np.mean([MSD_n_z_avg_1,MSD_n_z_avg_2,MSD_n_z_avg_3,MSD_n_z_avg_4,MSD_n_z_avg_5], axis=0)

    msd_n_xy_avg = np.mean([MSD_n_xy_avg_1,MSD_n_xy_avg_2,MSD_n_xy_avg_3,MSD_n_xy_avg_4,MSD_n_xy_avg_5], axis=0)
    msd_n_xz_avg = np.mean([MSD_n_xz_avg_1,MSD_n_xz_avg_2,MSD_n_xz_avg_3,MSD_n_xz_avg_4,MSD_n_xz_avg_5], axis=0)
    msd_n_yz_avg = np.mean([MSD_n_yz_avg_1,MSD_n_yz_avg_2,MSD_n_yz_avg_3,MSD_n_yz_avg_4,MSD_n_yz_avg_5], axis=0)

    msd_n_tot = msd_n_x_avg + msd_n_y_avg + msd_n_z_avg
    time_avg = time[0]

    import matplotlib.pyplot as plt


    # Convert the unit of the data 
    t = time_avg*1e-12         # Ps to sec 
    msd_x = msd_n_x_avg*1e-20     # A^2 to m^2
    msd_y = msd_n_y_avg*1e-20     # A^2 to m^2
    msd_z = msd_n_z_avg*1e-20     # A^2 to m^2
    msd_xy = msd_n_xy_avg*1e-20     # A^2 to m^2
    msd_yz = msd_n_yz_avg*1e-20     # A^2 to m^2
    msd_xz = msd_n_xz_avg*1e-20     # A^2 to m^2
    msd_tot = msd_n_tot*1e-20     # A^2 to m^2

    def func(x, a, b):
        return a*x + b

    sl = 1

    x = t[::sl]
    y_x = msd_x[::sl]
    y_y = msd_y[::sl]
    y_z = msd_z[::sl]
    y_xy = msd_xy[::sl]
    y_yz = msd_yz[::sl]
    y_xz = msd_xz[::sl]
    y_tot = msd_tot[::sl]

    popt_x, pcov_x = curve_fit(func, x, y_x)
    popt_y, pcov_y = curve_fit(func, x, y_y)
    popt_z, pcov_z = curve_fit(func, x, y_z)
    popt_xy, pcov_xy = curve_fit(func, x, y_xy)
    popt_yz, pcov_yz = curve_fit(func, x, y_yz)
    popt_xz, pcov_xz = curve_fit(func, x, y_xz)
    popt_tot, pcov_tot = curve_fit(func, x, y_tot)

    import numpy as np




    m_x = 1  # Dimensions to define for diff
    m_y = 1
    m_z = 1
    m_xy = 2
    m_yz = 2
    m_xz = 2
    m_tot = 3
    n_x = 2*m_x
    n_y = 2*m_y
    n_z = 2*m_z
    n_xy = 2*m_xy
    n_yz = 2*m_yz
    n_xz = 2*m_xz
    n_tot = 2*m_tot

    D_x = popt_x[0]/n_x
    D_y = popt_y[0]/n_y
    D_z = popt_z[0]/n_z
    D_xy = popt_xy[0]/n_xy
    D_yz = popt_yz[0]/n_yz
    D_xz = popt_xz[0]/n_xz
    D_tot = popt_tot[0]/n_tot

    D_x_std = np.sqrt(np.diag(pcov_x))[0]/n_x
    D_y_std = np.sqrt(np.diag(pcov_y))[0]/n_y
    D_z_std = np.sqrt(np.diag(pcov_z))[0]/n_z
    D_xy_std = np.sqrt(np.diag(pcov_xy))[0]/n_xy
    D_yz_std = np.sqrt(np.diag(pcov_yz))[0]/n_yz
    D_xz_std = np.sqrt(np.diag(pcov_xz))[0]/n_xz
    D_tot_std = np.sqrt(np.diag(pcov_tot))[0]/n_tot

    print("Diffusion coefficient x = %.3e [M^2/sec]" %(D_x))
    print("Diffusion coefficient y = %.3e [M^2/sec]" %(D_y))
    print("Diffusion coefficient z = %.3e [M^2/sec]" %(D_z))
    print("Diffusion coefficient tot = %.3e [M^2/sec]" %(D_tot))
    print("Stndard deviation x = %.3e [M^2/sec]" %(D_x_std))
    print("Stndard deviation y = %.3e [M^2/sec]" %(D_y_std))
    print("Stndard deviation z = %.3e [M^2/sec]" %(D_z_std))
    print("Stndard deviation tot = %.3e [M^2/sec]" %(D_tot_std))

    plt.rcParams["figure.facecolor"] = "w"
    plt.rcParams["figure.figsize"] = (7,7)
    plt.rcParams.update({'font.size': 14})
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.xlim(xlimit)
    
    if ylimt == True:
        plt.ylim(y_limit)
    
    plt.tick_params(bottom=True, top=True, left=True, right=True)
    plt.rcParams["axes.linewidth"] = 3
    plt.rcParams["patch.linewidth"] = 3
    plt.tick_params(axis="y",direction="in")
    plt.tick_params(axis="x",direction="in")
    plt.tick_params(width=3, length=4.5)
    plt.plot(x*1e+9, func(x, *popt_x)*(1e+10)**2,'--',color='black',linewidth=2)   #,label='$D_{tr}^{%s}$ = %.3e ($M^2$/sec)'%(temp,d))
    plt.plot(x*1e+9, func(x, *popt_y)*(1e+10)**2,'--',color='black',linewidth=2)   #,label='$D_{tr}^{%s}$ = %.3e ($M^2$/sec)'%(temp,d))
    plt.plot(x*1e+9, func(x, *popt_z)*(1e+10)**2,'--',color='black',linewidth=2)   #,label='$D_{tr}^{%s}$ = %.3e ($M^2$/sec)'%(temp,d))
    plt.plot(x*1e+9, func(x, *popt_tot)*(1e+10)**2,'--',color='black',linewidth=2) #,label='$D_{tr}^{%s}$ = %.3e ($M^2$/sec)'%(temp,d))

    plt.plot(x*1e+9,y_x*(1e+10)**2,label='$D_{X}^{*}$ = %.3e ($m^2$/s)'%(D_x), color='royalblue', linewidth=3, alpha=0.75)     #, label='Average MSD')
    plt.plot(x*1e+9,y_y*(1e+10)**2,label='$D_{Y}^{*}$ = %.3e ($m^2$/s)'%(D_y), color='darkorange',linewidth=3, alpha=0.75)     #, label='Average MSD')
    plt.plot(x*1e+9,y_z*(1e+10)**2,label='$D_{Z}^{*}$ = %.3e ($m^2$/s)'%(D_z), color='green',linewidth=3, alpha=0.75)     #, label='Average MSD')
    plt.plot(x*1e+9,y_tot*(1e+10)**2,label='$D_{Tot}^{*}$ = %.3e ($m^2$/s)'%(D_tot), color='red',linewidth=3, alpha=0.75) #, label='Average MSD')
    plt.xlabel(r"Time (ns)", fontsize=16, weight="bold")
    plt.ylabel(r"$MSD^{COM}$ ($\AA^2$)",fontsize=16)
    plt.legend(fancybox=False, framealpha=1, shadow=False, borderpad=1, frameon=True, edgecolor="black",fontsize="12")#, bbox_to_anchor=(1.60, 1.019))
    plt.savefig(output_file, bbox_inches='tight', dpi=600, transparent=False)
    plt.show()
    
    return D_x, D_y, D_z, D_xy, D_yz, D_xz, D_tot, D_x_std, D_y_std, D_z_std, D_xy_std, D_yz_std, D_xz_std, D_tot_std


def compute_msd_tracer(filee_list=['msd_na_nasi_523_1_0_5000.txt', 'msd_na_nasi_523_1_0_5000.txt', 'msd_na_nasi_523_1_0_5000.txt', 'msd_na_nasi_523_1_0_5000.txt', 'msd_na_nasi_523_1_0_5000.txt'], output_file='msd_glass_523k_avg_sigma.png', xlimit=[0, 2], ylimt=True, y_limit=[0, 900]):
    
    import numpy as np

    MSD_n_x_1 = []
    MSD_n_x_2 = []
    MSD_n_x_3 = []
    MSD_n_x_4 = []
    MSD_n_x_5 = []
    MSD_n_y_1 = []
    MSD_n_y_2 = []
    MSD_n_y_3 = []
    MSD_n_y_4 = []
    MSD_n_y_5 = []
    MSD_n_z_1 = []
    MSD_n_z_2 = []
    MSD_n_z_3 = []
    MSD_n_z_4 = []
    MSD_n_z_5 = []
    MSD_n_xy_1 = []
    MSD_n_xy_2 = []
    MSD_n_xy_3 = []
    MSD_n_xy_4 = []
    MSD_n_xy_5 = []
    MSD_n_xz_1 = []
    MSD_n_xz_2 = []
    MSD_n_xz_3 = []
    MSD_n_xz_4 = []
    MSD_n_xz_5 = []
    MSD_n_yz_1 = []
    MSD_n_yz_2 = []
    MSD_n_yz_3 = []
    MSD_n_yz_4 = []
    MSD_n_yz_5 = []
    time = []

    for i in range(1):
        if i == 0:

            file_data_1 = np.loadtxt(filee_list[0])
            file_data_2 = np.loadtxt(filee_list[1])
            file_data_3 = np.loadtxt(filee_list[2])
            file_data_4 = np.loadtxt(filee_list[3])
            file_data_5 = np.loadtxt(filee_list[4])
            
            t = file_data_1[:,0]
            time.append(t)
        
            netMSD_x_1 = file_data_1[:,1]
            netMSD_x_2 = file_data_2[:,1]
            netMSD_x_3 = file_data_3[:,1]
            netMSD_x_4 = file_data_4[:,1]
            netMSD_x_5 = file_data_5[:,1]

            netMSD_y_1 = file_data_1[:,2]
            netMSD_y_2 = file_data_2[:,2]
            netMSD_y_3 = file_data_3[:,2]
            netMSD_y_4 = file_data_4[:,2]
            netMSD_y_5 = file_data_5[:,2]

            netMSD_z_1 = file_data_1[:,3]
            netMSD_z_2 = file_data_2[:,3]   
            netMSD_z_3 = file_data_3[:,3]
            netMSD_z_4 = file_data_4[:,3]
            netMSD_z_5 = file_data_5[:,3]

            netMSD_xy_1 = file_data_1[:,4]
            netMSD_xy_2 = file_data_2[:,4]
            netMSD_xy_3 = file_data_3[:,4]
            netMSD_xy_4 = file_data_4[:,4]
            netMSD_xy_5 = file_data_5[:,4]

            netMSD_xz_1 = file_data_1[:,5]
            netMSD_xz_2 = file_data_2[:,5]
            netMSD_xz_3 = file_data_3[:,5]
            netMSD_xz_4 = file_data_4[:,5]
            netMSD_xz_5 = file_data_5[:,5]
            
            netMSD_yz_1 = file_data_1[:,5]
            netMSD_yz_2 = file_data_2[:,5]
            netMSD_yz_3 = file_data_3[:,5]
            netMSD_yz_4 = file_data_4[:,5]
            netMSD_yz_5 = file_data_5[:,5]

            MSD_n_x_1.append(netMSD_x_1)
            MSD_n_x_2.append(netMSD_x_2)
            MSD_n_x_3.append(netMSD_x_3)
            MSD_n_x_4.append(netMSD_x_4)
            MSD_n_x_5.append(netMSD_x_5)

            MSD_n_y_1.append(netMSD_y_1)
            MSD_n_y_2.append(netMSD_y_2)
            MSD_n_y_3.append(netMSD_y_3)
            MSD_n_y_4.append(netMSD_y_4)
            MSD_n_y_5.append(netMSD_y_5)

            MSD_n_z_1.append(netMSD_z_1)
            MSD_n_z_2.append(netMSD_z_2)
            MSD_n_z_3.append(netMSD_z_3)
            MSD_n_z_4.append(netMSD_z_4)
            MSD_n_z_5.append(netMSD_z_5)

            MSD_n_xy_1.append(netMSD_xy_1)
            MSD_n_xy_2.append(netMSD_xy_2)
            MSD_n_xy_3.append(netMSD_xy_3)
            MSD_n_xy_4.append(netMSD_xy_4)
            MSD_n_xy_5.append(netMSD_xy_5)

            MSD_n_xz_1.append(netMSD_xz_1)
            MSD_n_xz_2.append(netMSD_xz_2)
            MSD_n_xz_3.append(netMSD_xz_3)
            MSD_n_xz_4.append(netMSD_xz_4)
            MSD_n_xz_5.append(netMSD_xz_5)

            MSD_n_yz_1.append(netMSD_yz_1)
            MSD_n_yz_2.append(netMSD_yz_2)
            MSD_n_yz_3.append(netMSD_yz_3)
            MSD_n_yz_4.append(netMSD_yz_4)
            MSD_n_yz_5.append(netMSD_yz_5)

    MSD_n_x_avg_1 = np.mean(MSD_n_x_1,axis=0)
    MSD_n_x_avg_2 = np.mean(MSD_n_x_2,axis=0)
    MSD_n_x_avg_3 = np.mean(MSD_n_x_3,axis=0)
    MSD_n_x_avg_4 = np.mean(MSD_n_x_4,axis=0)
    MSD_n_x_avg_5 = np.mean(MSD_n_x_5,axis=0)

    MSD_n_y_avg_1 = np.mean(MSD_n_y_1,axis=0)
    MSD_n_y_avg_2 = np.mean(MSD_n_y_2,axis=0)
    MSD_n_y_avg_3 = np.mean(MSD_n_y_3,axis=0)
    MSD_n_y_avg_4 = np.mean(MSD_n_y_4,axis=0)
    MSD_n_y_avg_5 = np.mean(MSD_n_y_5,axis=0)

    MSD_n_z_avg_1 = np.mean(MSD_n_z_1,axis=0)
    MSD_n_z_avg_2 = np.mean(MSD_n_z_2,axis=0)
    MSD_n_z_avg_3 = np.mean(MSD_n_z_3,axis=0)
    MSD_n_z_avg_4 = np.mean(MSD_n_z_4,axis=0)
    MSD_n_z_avg_5 = np.mean(MSD_n_z_5,axis=0)

    MSD_n_xy_avg_1 = np.mean(MSD_n_xy_1,axis=0)
    MSD_n_xy_avg_2 = np.mean(MSD_n_xy_2,axis=0)
    MSD_n_xy_avg_3 = np.mean(MSD_n_xy_3,axis=0)
    MSD_n_xy_avg_4 = np.mean(MSD_n_xy_4,axis=0)
    MSD_n_xy_avg_5 = np.mean(MSD_n_xy_5,axis=0)

    MSD_n_xz_avg_1 = np.mean(MSD_n_xz_1,axis=0)
    MSD_n_xz_avg_2 = np.mean(MSD_n_xz_2,axis=0)
    MSD_n_xz_avg_3 = np.mean(MSD_n_xz_3,axis=0)
    MSD_n_xz_avg_4 = np.mean(MSD_n_xz_4,axis=0)
    MSD_n_xz_avg_5 = np.mean(MSD_n_xz_5,axis=0)

    MSD_n_yz_avg_1 = np.mean(MSD_n_yz_1,axis=0)
    MSD_n_yz_avg_2 = np.mean(MSD_n_yz_2,axis=0)
    MSD_n_yz_avg_3 = np.mean(MSD_n_yz_3,axis=0)
    MSD_n_yz_avg_4 = np.mean(MSD_n_yz_4,axis=0)
    MSD_n_yz_avg_5 = np.mean(MSD_n_yz_5,axis=0)

    msd_n_x_avg = np.mean([MSD_n_x_avg_1,MSD_n_x_avg_2,MSD_n_x_avg_3,MSD_n_x_avg_4,MSD_n_x_avg_5], axis=0)
    msd_n_y_avg = np.mean([MSD_n_y_avg_1,MSD_n_y_avg_2,MSD_n_y_avg_3,MSD_n_y_avg_4,MSD_n_y_avg_5], axis=0)
    msd_n_z_avg = np.mean([MSD_n_z_avg_1,MSD_n_z_avg_2,MSD_n_z_avg_3,MSD_n_z_avg_4,MSD_n_z_avg_5], axis=0)

    msd_n_xy_avg = np.mean([MSD_n_xy_avg_1,MSD_n_xy_avg_2,MSD_n_xy_avg_3,MSD_n_xy_avg_4,MSD_n_xy_avg_5], axis=0)
    msd_n_xz_avg = np.mean([MSD_n_xz_avg_1,MSD_n_xz_avg_2,MSD_n_xz_avg_3,MSD_n_xz_avg_4,MSD_n_xz_avg_5], axis=0)
    msd_n_yz_avg = np.mean([MSD_n_yz_avg_1,MSD_n_yz_avg_2,MSD_n_yz_avg_3,MSD_n_yz_avg_4,MSD_n_yz_avg_5], axis=0)

    msd_n_tot = msd_n_x_avg + msd_n_y_avg + msd_n_z_avg
    time_avg = time[0]

    import matplotlib.pyplot as plt


    # Convert the unit of the data 
    t = time_avg*1e-12         # Ps to sec 
    msd_x = msd_n_x_avg*1e-20     # A^2 to m^2
    msd_y = msd_n_y_avg*1e-20     # A^2 to m^2
    msd_z = msd_n_z_avg*1e-20     # A^2 to m^2
    msd_xy = msd_n_xy_avg*1e-20     # A^2 to m^2
    msd_yz = msd_n_yz_avg*1e-20     # A^2 to m^2
    msd_xz = msd_n_xz_avg*1e-20     # A^2 to m^2
    msd_tot = msd_n_tot*1e-20     # A^2 to m^2

    def func(x, a, b):
        return a*x + b

    sl = 1

    x = t[::sl]
    y_x = msd_x[::sl]
    y_y = msd_y[::sl]
    y_z = msd_z[::sl]
    y_xy = msd_xy[::sl]
    y_yz = msd_yz[::sl]
    y_xz = msd_xz[::sl]
    y_tot = msd_tot[::sl]

    popt_x, pcov_x = curve_fit(func, x, y_x)
    popt_y, pcov_y = curve_fit(func, x, y_y)
    popt_z, pcov_z = curve_fit(func, x, y_z)
    popt_xy, pcov_xy = curve_fit(func, x, y_xy)
    popt_yz, pcov_yz = curve_fit(func, x, y_yz)
    popt_xz, pcov_xz = curve_fit(func, x, y_xz)
    popt_tot, pcov_tot = curve_fit(func, x, y_tot)

    import numpy as np

    m_x = 1  # Dimensions to define for diff
    m_y = 1
    m_z = 1
    m_xy = 2
    m_yz = 2
    m_xz = 2
    m_tot = 3
    n_x = 2*m_x
    n_y = 2*m_y
    n_z = 2*m_z
    n_xy = 2*m_xy
    n_yz = 2*m_yz
    n_xz = 2*m_xz
    n_tot = 2*m_tot

    D_x = popt_x[0]/n_x
    D_y = popt_y[0]/n_y
    D_z = popt_z[0]/n_z
    D_xy = popt_xy[0]/n_xy
    D_yz = popt_yz[0]/n_yz
    D_xz = popt_xz[0]/n_xz
    D_tot = popt_tot[0]/n_tot

    D_x_std = np.sqrt(np.diag(pcov_x))[0]/n_x
    D_y_std = np.sqrt(np.diag(pcov_y))[0]/n_y
    D_z_std = np.sqrt(np.diag(pcov_z))[0]/n_z
    D_xy_std = np.sqrt(np.diag(pcov_xy))[0]/n_xy
    D_yz_std = np.sqrt(np.diag(pcov_yz))[0]/n_yz
    D_xz_std = np.sqrt(np.diag(pcov_xz))[0]/n_xz
    D_tot_std = np.sqrt(np.diag(pcov_tot))[0]/n_tot

    print("Diffusion coefficient x = %.3e [M^2/sec]" %(D_x))
    print("Diffusion coefficient y = %.3e [M^2/sec]" %(D_y))
    print("Diffusion coefficient z = %.3e [M^2/sec]" %(D_z))
    print("Diffusion coefficient tot = %.3e [M^2/sec]" %(D_tot))
    print("Stndard deviation x = %.3e [M^2/sec]" %(D_x_std))
    print("Stndard deviation y = %.3e [M^2/sec]" %(D_y_std))
    print("Stndard deviation z = %.3e [M^2/sec]" %(D_z_std))
    print("Stndard deviation tot = %.3e [M^2/sec]" %(D_tot_std))

    plt.rcParams["figure.facecolor"] = "w"
    plt.rcParams["figure.figsize"] = (7,7)
    plt.rcParams.update({'font.size': 14})
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.xlim(xlimit)
    
    if ylimt == True:
        plt.ylim(y_limit)
    
    plt.tick_params(bottom=True, top=True, left=True, right=True)
    plt.rcParams["axes.linewidth"] = 3
    plt.rcParams["patch.linewidth"] = 3
    plt.tick_params(axis="y",direction="in")
    plt.tick_params(axis="x",direction="in")
    plt.tick_params(width=3, length=4.5)
    plt.plot(x*1e+9, func(x, *popt_x)*(1e+10)**2,'--',color='black',linewidth=2)   #,label='$D_{tr}^{%s}$ = %.3e ($M^2$/sec)'%(temp,d))
    plt.plot(x*1e+9, func(x, *popt_y)*(1e+10)**2,'--',color='black',linewidth=2)   #,label='$D_{tr}^{%s}$ = %.3e ($M^2$/sec)'%(temp,d))
    plt.plot(x*1e+9, func(x, *popt_z)*(1e+10)**2,'--',color='black',linewidth=2)   #,label='$D_{tr}^{%s}$ = %.3e ($M^2$/sec)'%(temp,d))
    plt.plot(x*1e+9, func(x, *popt_tot)*(1e+10)**2,'--',color='black',linewidth=2) #,label='$D_{tr}^{%s}$ = %.3e ($M^2$/sec)'%(temp,d))

    plt.plot(x*1e+9,y_x*(1e+10)**2,label='$D_{X}^{*}$ = %.3e ($m^2$/s)'%(D_x), color='royalblue', linewidth=3, alpha=0.75)     #, label='Average MSD')
    plt.plot(x*1e+9,y_y*(1e+10)**2,label='$D_{Y}^{*}$ = %.3e ($m^2$/s)'%(D_y), color='darkorange',linewidth=3, alpha=0.75)     #, label='Average MSD')
    plt.plot(x*1e+9,y_z*(1e+10)**2,label='$D_{Z}^{*}$ = %.3e ($m^2$/s)'%(D_z), color='green',linewidth=3, alpha=0.75)     #, label='Average MSD')
    plt.plot(x*1e+9,y_tot*(1e+10)**2,label='$D_{Tot}^{*}$ = %.3e ($m^2$/s)'%(D_tot), color='red',linewidth=3, alpha=0.75) #, label='Average MSD')
    plt.xlabel(r"Time (ns)", fontsize=16, weight="bold")
    plt.ylabel(r"$MSD^{COM}$ ($\AA^2$)",fontsize=16)
    plt.legend(fancybox=False, framealpha=1, shadow=False, borderpad=1, frameon=True, edgecolor="black",fontsize="12")#, bbox_to_anchor=(1.60, 1.019))
    plt.savefig(output_file, bbox_inches='tight', dpi=600, transparent=False)
    plt.show()
    
    return D_x, D_y, D_z, D_xy, D_yz, D_xz, D_tot, D_x_std, D_y_std, D_z_std, D_xy_std, D_yz_std, D_xz_std, D_tot_std