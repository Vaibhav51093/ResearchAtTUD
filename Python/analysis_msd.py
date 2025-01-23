#-------------------------------------------------------------------------------------------------------------------------------#
# Official work of Vaibhav Arun Deshmukh,                                                                                       #
# Module to calculate the following (bulk, gb and interface)                                                                    #
#       1. Tracer diffusion                                                                                                     #
#       2. Ein. diffusion                                                                                                       #   
#       3. Ionic conductivity                                                                                                   # 
#       4. Resistance                                                                                                           #
#-------------------------------------------------------------------------------------------------------------------------------#  

from fileinput import filename
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
#from scipy import stats as st
import tarfile
import lammps_logfile
import os
import codecs
utf8reader = codecs.getreader('utf-8')
import os
import numpy as np 
from scipy import integrate
from scipy.constants import codata
import pandas as pd
from ase.io import read, write
from scipy.optimize import curve_fit
#from pyiron import Project
#from pyiron import ase_to_pyiron, pyiron_to_ase
import shutil
import glob
from scipy import stats
from scipy.optimize import curve_fit
import numpy as np 
import scipy.constants as const

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev


#--------------------------------------------------------MSD-module-------------------------------------------------------------#

class msd():

    def dff_raw(file="cryt_msd_523k_5.txt",dim=3,temp=523,y_val=1):

        from scipy.optimize import curve_fit
        import numpy as np 
        # Import data and select as per file 
        file_data = np.loadtxt(file)
        x = file_data[:,0]                 # Time 
        y_avg = (file_data[:,y_val])#**2)/162 # 17 to 20 avg, x, y, z, 3, 1,1,1 dimensions
        m = dim  # Dimensions to define for diff
        n = 2*m 

        # Smoothening the data 
        xy = np.column_stack((x, y_avg))
        z = pd.DataFrame(xy).groupby(0, as_index=False)[1].mean().values

        x_s = z[:, 0] 
        y_s = z[:, 1]
    
        t = x_s*1e-12       # to sec 
        msd = y_s*1e-20     # A^2 to m^2

        # Define the functiion to fit the curve 
        def func(x, a, b):
            return a*x + b

        x_1 = t
        y = msd
        
        popt, pcov = curve_fit(func, x_1, y)
        #print(popt[0]/2)  
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_1,y)

        print("Diffusion coefficient = %.3e [M^2/sec]" %(slope/n))

        d = slope/n

        return 

    def dff_coefficient(file="cryt_msd_523k_5.txt",dim=3,temp=523,save=False,file_2="avg_msd_773k_5.png",y_val=1,sl=1):
        
        from scipy.optimize import curve_fit
        import numpy as np 
        # Import data and select as per file 
        file_data = np.loadtxt(file)
        x = file_data[:,0]                 # Time 
        y_avg = (file_data[:,y_val])#**2)/162 # 17 to 20 avg, x, y, z, 3, 1,1,1 dimensions
        m = dim  # Dimensions to define for diff
        n = 2*m 

        # Smoothening the data 
        xy = np.column_stack((x, y_avg))
        z = pd.DataFrame(xy).groupby(0, as_index=False)[1].mean().values

        x_s = z[:, 0] 
        y_s = z[:, 1]
    
        t = x_s*1e-12       # to sec 
        msd = y_s*1e-20     # A^2 to m^2

        # Define the functiion to fit the curve 
        def func(x, a, b):
            return a*x + b

        x_1 = t[::sl]
        y = msd[::sl]
        
        popt, pcov = curve_fit(func, x_1, y)
        #print(popt[0]/2)  
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_1,y)

        print("Diffusion coefficient = %.3e [M^2/sec]" %(slope/n))

        d = slope/n

        
        import matplotlib.pyplot as plt
        plt.rcParams["figure.facecolor"] = "w"
        plt.rcParams["figure.figsize"] = (7,7)
        plt.rcParams.update({'font.size': 12})
        plt.plot(x_1, func(x_1, *popt),'--')#,label='$D_{tr}^{%s}$ = %.3e ($M^2$/sec)'%(temp,d))
        plt.plot(x_1,y,label='$D_{tr}^{%s}$ = %.3e ($M^2$/sec)'%(temp,d))#, label='Average MSD')
        #plt.xlabel("MD steps")
        #plt.ylabel("MSD")
        plt.grid(which='both')
        plt.grid(True, linestyle='--')
        #plt.title("Tracer Diffusion Coefficient")
        plt.xlabel(r"Time ($Sec$)", fontsize=14)
        plt.ylabel("$MSD\ (M^2)$",fontsize=14)
        plt.legend(fancybox=False, framealpha=1, shadow=False, borderpad=1, frameon=True, edgecolor="black")#, bbox_to_anchor=(1.60, 1.019))
        if save == True:
            plt.savefig(file_2, bbox_inches='tight', dpi=600, transparent=False)
        plt.show 

        return d, x, y

    def roll_mean(data_1, data_2, window):
        avg_data_1 = []
        avg_data_2 = []
        for i,j in zip(range(len(data_1)-window+1),range(len(data_2)-window+1)):
            avg_data_1.append(np.mean(data_1[i:i+window]))
            avg_data_2.append(np.mean(data_2[j:j+window]))
        return avg_data_1,avg_data_2

    def roll_mean_1(data_1, window):
        avg_data_1 = []
        for i in zip(range(len(data_1)-window+1)):
            avg_data_1.append(np.mean(data_1[i:i+window]))
        return avg_data_1


    def smooth_msd_data(x, y):
        """
        Smooths the data from the smoothed msd function. The data consists
        of multiple msd runs but the data is unordered. This function will
        order the data and average all y values with equivalent x values.
        Args:
            x (:py:attr:`array_like`): Time data.
            y (:py:attr:`array_like`): MSD data.
        
        Returns:
            z (:py:attr:`array_like`): Time / MSD data.
        """
        xy = np.column_stack((x, y))
        z = pd.DataFrame(xy).groupby(0, as_index=False)[1].mean().values
        return z[:, 0], z[:, 1]

    def fit_arrhenius(temps, diffusivities):
        """
        Returns Ea, c, standard error of Ea from the Arrhenius fit:
        D = c * exp(-Ea/kT)
        Args:
            temps ([float]): A sequence of temperatures. units: K
            diffusivities ([float]): A sequence of diffusivities (e.g.,
                from DiffusionAnalyzer.diffusivity). units: cm^2/s
        """
        import scipy.constants as const
        from scipy.stats import linregress
        t_1 = 1 / np.array(temps)
        diff = [] 
        for i in diffusivities:
            diff.append(i*10000)
        logd = np.log(diff)

        # Do a least squares regression of log(D) vs 1/T
        a = np.array([t_1, np.ones(len(temps))]).T
        w, res, _, _ = np.linalg.lstsq(a, logd, rcond=None)
        w = np.array(w)
        n = len(temps)
        
        if n > 2:
            std_Ea = (res[0] / (n - 2) / (n * np.var(t_1))) ** 0.5 * const.k / const.e
        else:
            std_Ea = None
        return -w[0] * const.k / const.e, np.exp(w[1]), std_Ea

    def get_extrapolated_diffusivity(temps, diffusivities, new_temp):
        """
            Returns (Arrhenius) extrapolated diffusivity at new_temp
        Args:
            temps ([float]): A sequence of temperatures. units: K
            diffusivities ([float]): A sequence of diffusivities (e.g.,
            from DiffusionAnalyzer.diffusivity). units: cm^2/s
            new_temp (float): desired temperature. units: K
        Returns:
            (float) Diffusivity at extrapolated temp in cm^2/s.
        """
        import scipy.constants as const
        Ea, c, _ = msd.fit_arrhenius(temps, diffusivities)
        return (c * np.exp(-Ea / (const.k / const.e * new_temp)))*0.0001

    def conductivity(charge=1, no_na_atoms=1020,diffu=5.510961250898897e-13,volume=120, temp=303):

        x_co = (1000 * no_na_atoms / ((volume * (1e-24))) * (charge)**2 * (const.e) ** 2 / (const.k * temp)) * (diffu) * 10000 # Return cond in mS/cm

        return x_co

    def calculate_ionic_conductivity(charge, no_of_ions, diffusion_coefficient, volume, temperature):
        # Convert volume from metal units to cubic meters
        volume_m3 = volume * 1e-30

        # Calculate the molar conductivity
        molar_conductivity = (charge ** 2 * const.e ** 2 * no_of_ions * diffusion_coefficient) / (const.k * temperature)

        # Convert molar conductivity to ionic conductivity (in mS/cm)
        ionic_conductivity = molar_conductivity * 1000 / (const.N_A * volume_m3) * 1e2

        return ionic_conductivity  #mS/cm


    def get_extrapolated_conductivity(temps, diffusivities, new_temp,charge=1, no_na_atoms=1020,volume=120):

        x_con = msd.get_extrapolated_diffusivity(temps, diffusivities, new_temp)

        x_co = (1000 * no_na_atoms / ((volume * (1e-24))) * (charge)**2 * (const.e) ** 2 / (const.k * new_temp)) * (x_con) * 10000 # mS/cm 

        return x_co 

    def get_arrhenius_plot(temps, diffusivities, diffusivity_errors=None,save=True,file="activation_energy_cryst_amorph_latest.png", **kwargs):
        r"""
        Returns an Arrhenius plot.
        Args:
            temps ([float]): A sequence of temperatures.
            diffusivities ([float]): A sequence of diffusivities (e.g.,
            from DiffusionAnalyzer.diffusivity).
            diffusivity_errors ([float]): A sequence of errors for the
            diffusivities. If None, no error bar is plotted.
            \\*\\*kwargs:
                Any keyword args supported by matplotlib.pyplot.plot.
        Returns:
            A matplotlib.pyplot object. Do plt.show() to show the plot.
        """
        Ea, c, _ = msd.fit_arrhenius(temps, diffusivities)     # Activation energy 

        # log10 of the arrhenius fit
        arr = (c * np.exp(-Ea / (const.k / const.e * np.array(temps))))*0.0001 # get log of diffusion
        #diff_2 = []   # m2/sec
        #for i in arr: 
        #    diff_2.append(i/10000)
        
        std_arr = diffusivity_errors
        #for i in arr:
        #    std_arr.append(np.std(i))

        t_1 = 1000 / np.array(temps)           # temperature
        plusminus = u"\u00B1"
        import matplotlib.pyplot as plt

        plt.rcParams["figure.facecolor"] = "w"
        plt.rcParams["figure.figsize"] = (7,7)
        plt.rcParams.update({'font.size': 14})
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["axes.labelweight"] = "bold"
        plt.rcParams["figure.figsize"] = (7,7)
        plt.tick_params(bottom=True, top=True, left=True, right=True)
        plt.tick_params(axis="y",direction="in")
        plt.tick_params(axis="x",direction="in")
        plt.rcParams["axes.linewidth"] = 3
        plt.rcParams["patch.linewidth"] = 3
        plt.tick_params(width=3, length=4.5)
        plt.errorbar(t_1,arr,yerr=std_arr,linestyle='--',fmt = 'o',color = 'blue', ecolor = 'blue', elinewidth = 3, capsize=5,label='$E_a$ = {:0.4f} {} {:0.3f} eV'.format(Ea,f'{plusminus}',_))
        #plt.plot(t_1, np.exp(result.slope * np.asarray(t_1) + result.intercept), '--', color = color)
        plt.yscale("log")
        #plt.title("Arrhenius Plot Crystal vs Amorphous (Ordered)")
        plt.xlabel("1000/$Temperature \ [1/K]$", fontsize=14)
        plt.ylabel("$log(D\ (m^2/s)$)", fontsize=14)
        plt.legend(fancybox=False, framealpha=1, shadow=False, borderpad=1, frameon=True, edgecolor="black", loc='upper right')#, bbox_to_anchor=(1.60, 1.019))
        if save == True:
            plt.savefig(file, bbox_inches='tight', dpi=600, transparent=False)
        plt.show()


    def activation_energy(temp = [523,573,623,673,723,773], d_523 = [],d_573 = [],d_623 = [], d_673 = [], d_723 = [], d_773 = [], color = 'blue', new=True, new_temp=300):#,extrapolate=False, new_temp=300):
        
        temp_new = []    
        t_1 = []                   # 1/K
        for i in temp:
            temp_new.append((1/i))
            t_1.append((1/i)*1000)
        d_cryt = [np.mean(d_523),np.mean(d_573),np.mean(d_623),np.mean(d_673),np.mean(d_723),np.mean(d_773)]
        d_cryt_err = [np.std(d_523),np.std(d_523),np.std(d_623),np.std(d_673),np.std(d_723),np.std(d_773)]
        
        from scipy.stats import linregress
        from scipy.constants import R
        result = linregress(temp_new,np.log(d_cryt))
        activation_energy = - R * result.slope
        preexponential_factor = np.exp(result.intercept)
        print(f"for crystal; A = {preexponential_factor}; E_a = {activation_energy*10**(-3)} kJ/mol; E_a = {0.01*activation_energy*10**(-3)} eV")
        
        plusminus = u"\u00B1"
        #print(f'{plusminus}')

        err = (0.01)*result.stderr*R*10**(-3)
         
        import matplotlib.pyplot as plt
        plt.rcParams["figure.facecolor"] = "w"
        plt.rcParams["figure.figsize"] = (7,7)
        plt.rcParams.update({'font.size': 12})
        plt.tick_params(bottom=True, top=True, left=True, right=True)
        plt.tick_params(axis="y",direction="in")
        plt.tick_params(axis="x",direction="in")
        plt.errorbar(t_1,d_cryt,yerr=d_cryt_err,fmt = 'o',color = 'blue', ecolor = 'blue', elinewidth = 2, capsize=5,label='$E_a$ = {:0.4f} {} {:0.3f} eV'.format(0.01*activation_energy*10**(-3),f'{plusminus}',err))
        #plt.scatter(temp_new, d_cryt,marker='o', facecolors='none', label='E_a = {:0.4f} eV (Crystal)'.format(0.01*activation_energy*10**(-3)), color = 'blue')
        plt.plot(t_1, np.exp(result.slope * np.asarray(temp_new) + result.intercept), '--', color = color)
        plt.yscale("log")
        #plt.title("Arrhenius Plot Crystal vs Amorphous (Ordered)")
        plt.xlabel("$Temperature \ [1/K]$", fontsize=14)
        plt.ylabel("$Diffusion\ [m^2/Sec]$", fontsize=14)
        plt.legend(fancybox=False, framealpha=1, shadow=False, borderpad=1, frameon=True, edgecolor="black", loc='upper right')#, bbox_to_anchor=(1.60, 1.019))
        plt.savefig("activation_energy_cryst_amorph_latest.png", bbox_inches='tight', dpi=600, transparent=False)
        plt.show()

        if new == True:#
            x = preexponential_factor * np.exp(- activation_energy / (const.R * new_temp))
            print('Extrapolated diffusivity = %.3e [m^2/Sec]'%(x))       

        return 0.01*activation_energy*10**(-3),x 


    def msd(file,file_2='msd_net.txt'):

        # Simulation Parameters
        POTIM = step = 1.00                 # Time step in fs
        NBLOCK = dump_write = 100           # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        # MSD Evaluation
        Stepsize = step_ovito = 10
        pipeline = import_file(file)
        # Print the list of input particle type

        Stepsize = step_ovito = 10
        pipeline = import_file(file)
        # Print the list of input particle type

        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)

        # Extract Atom Types
        AtomTypes        = []
        data = pipeline.compute()
        for Atom in data.particles['Particle Type']:
            if data.particles['Particle Type'].type_by_id(Atom).name not in AtomTypes:
                AtomTypes.append(data.particles['Particle Type'].type_by_id(Atom).name)


        # Extract number of Atoms for each Atom Type
        AtomTypesCount   = []
        for Type in AtomTypes:
            count = 0
            for Atom in data.particles['Particle Type']:
                if data.particles['Particle Type'].type_by_id(Atom).name == Type:
                    count += 1
            AtomTypesCount.append(count)


        # Set up Displacement Modifier and append to pipeline
        displ_mod    = CalculateDisplacementsModifier(
            use_frame_offset  = True,
            frame_offset      = -Stepsize)
        pipeline.modifiers.append(displ_mod)


        # Set up empty array to store intermediate displacements
        Displacements = []
        for Atom in data.particles['Particle Type']:
            Displacements.append([0.0,0.0,0.0])


        # Extract Displacements
        Time_and_MSD = [] #Structure will be:  [ Time, Li_Total, [Li_X,Li_Y,Li_Z], P_Total, [P_X,P_Y,P_Z], ... ]
        print("Starting MSD calculation for complete Trajectory")
        for frame in range(Stepsize,pipeline.source.num_frames,Stepsize):

            # Just to see a bit of update/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization
            if frame % 1000 == 0:
                print("    Computing frame %d of %d" % (frame,pipeline.source.num_frames),flush=True)


            # Set up temporary List that is added to Time_and_MSD list after each evaluation
            # Structure is:  [ Time, Li_Total, [Li_X,Li_Y,Li_Z], P_Total, [P_X,P_Y,P_Z], ... ]
            Time_and_MSD_temporary = [frame*NBLOCK*POTIM/1000]
            for Atom in AtomTypes:
                Time_and_MSD_temporary.append(0.0)
                Time_and_MSD_temporary.append([0.0,0.0,0.0])


            # Evaluate pipeline to let the modifier compute the Displacements of the current frame:
            data = pipeline.compute(frame)

            # Each iteration: Update the displacement vector of each atom, check its atom type and add the MSD to the Time_and_MSD_temporary list.
            for (AtomID, Atom, displ) in zip(range(0, sum(AtomTypesCount)), data.particles['Particle Type'], data.particles['Displacement']):
                Displacements[AtomID] = np.add(Displacements[AtomID],displ)
                for ListPosition,Type in enumerate(AtomTypes):
                    if data.particles['Particle Type'].type_by_id(Atom).name == Type:
                        Time_and_MSD_temporary[1+ListPosition*2]      += (Displacements[AtomID][0]**2 + Displacements[AtomID][1]**2 + Displacements[AtomID][2]**2) / AtomTypesCount[ListPosition]
                        Time_and_MSD_temporary[1+ListPosition*2+1][0] += (Displacements[AtomID][0]**2                                                            ) / AtomTypesCount[ListPosition]
                        Time_and_MSD_temporary[1+ListPosition*2+1][1] += (                              Displacements[AtomID][1]**2                              ) / AtomTypesCount[ListPosition]
                        Time_and_MSD_temporary[1+ListPosition*2+1][2] += (                                                            Displacements[AtomID][2]**2) / AtomTypesCount[ListPosition]
            Time_and_MSD.append(Time_and_MSD_temporary)


        # Header File
        header = '#MSD evaluated over all available data. Time in ps and MSDs in Ang^2\n#Time'
        for Type in AtomTypes:
            header = header + '\t  ' + Type + '\t  ' + Type + '_X\t  ' + Type + '_Y\t  ' + Type + '_Z'
        header = header + '\n  0.000\t'
        for Type in AtomTypes:
            header = header + '  0.000\t' + '  0.000\t' + '  0.000\t' + '  0.000\t'
        header = header + '\n'

        with open(file_2, 'w') as f:
            f.write(header)
            for DataPoints in Time_and_MSD:
                f.write('{:7.3f}'.format(DataPoints[0])) # Time
                for ListPosition in range(0,len(AtomTypes)):
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2  ]   ))
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2+1][0]))
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2+1][1]))
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2+1][2]))
                f.write('\n')
        #:6.2f

    def ein_msd(file,file_2='msd_net.txt'):
        
        # Simulation Parameters
        POTIM = step = 1.00                 # Time step in fs
        NBLOCK = dump_write = 100           # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        # MSD Evaluation
        Stepsize = step_ovito = 10
        pipeline = import_file(file)
        # Print the list of input particle type

        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)

        # Extract Atom Types
        AtomTypes        = []
        data = pipeline.compute()
        for Atom in data.particles['Particle Type']:
            if data.particles['Particle Type'].type_by_id(Atom).name not in AtomTypes:
                AtomTypes.append(data.particles['Particle Type'].type_by_id(Atom).name)


        # Extract number of Atoms for each Atom Type
        AtomTypesCount   = []
        for Type in AtomTypes:
            count = 0
            for Atom in data.particles['Particle Type']:
                if data.particles['Particle Type'].type_by_id(Atom).name == Type:
                    count += 1
            AtomTypesCount.append(count)


        # Set up Displacement Modifier and append to pipeline
        displ_mod    = CalculateDisplacementsModifier(
            use_frame_offset  = True,
            frame_offset      = -Stepsize)
        pipeline.modifiers.append(displ_mod)


        # Set up empty array to store intermediate displacements
        Displacements = []
        for Atom in data.particles['Particle Type']:
            Displacements.append([0.0,0.0,0.0])


        # Extract Displacements
        Time_and_MSD = [] #Structure will be:  [ Time, Li_Total, [Li_X,Li_Y,Li_Z], P_Total, [P_X,P_Y,P_Z], ... ]
        print("Starting MSD calculation for complete Trajectory")
        for frame in range(Stepsize,pipeline.source.num_frames,Stepsize):

            # Just to see a bit of update
            if frame % 1000 == 0:
                print("    Computing frame %d of %d" % (frame,pipeline.source.num_frames),flush=True)


            # Set up temporary List that is added to Time_and_MSD list after each evaluation
            # Structure is:  [ Time, Li_Total, [Li_X,Li_Y,Li_Z], P_Total, [P_X,P_Y,P_Z], ... ]
            Time_and_MSD_temporary = [frame*NBLOCK*POTIM/1000]
            for Atom in AtomTypes:
                Time_and_MSD_temporary.append(0.0)
                Time_and_MSD_temporary.append([0.0,0.0,0.0])


            # Evaluate pipeline to let the modifier compute the Displacements of the current frame:
            data = pipeline.compute(frame)

            # Each iteration: Update the displacement vector of each atom, check its atom type and add the MSD to the Time_and_MSD_temporary list.
            for (AtomID, Atom, displ) in zip(range(0, sum(AtomTypesCount)), data.particles['Particle Type'], data.particles['Displacement']):
                Displacements[AtomID] = np.add(Displacements[AtomID],displ)
                for ListPosition,Type in enumerate(AtomTypes):
                    if data.particles['Particle Type'].type_by_id(Atom).name == Type:
                        Time_and_MSD_temporary[1+ListPosition*2]      += (Displacements[AtomID][0] + Displacements[AtomID][1] + Displacements[AtomID][2])       # / AtomTypesCount[ListPosition]
                        Time_and_MSD_temporary[1+ListPosition*2+1][0] += (Displacements[AtomID][0]                                                            ) # / AtomTypesCount[ListPosition]
                        Time_and_MSD_temporary[1+ListPosition*2+1][1] += (                              Displacements[AtomID][1]                              ) # / AtomTypesCount[ListPosition]
                        Time_and_MSD_temporary[1+ListPosition*2+1][2] += (                                                            Displacements[AtomID][2]) # / AtomTypesCount[ListPosition]
            Time_and_MSD.append(Time_and_MSD_temporary)


        # Header File
        header = '#MSD evaluated over all available data. Time in ps and MSDs in Ang^2\n#Time'
        for Type in AtomTypes:
            header = header + '\t  ' + Type + '\t  ' + Type + '_X\t  ' + Type + '_Y\t  ' + Type + '_Z'
        header = header + '\n  0.000\t'
        for Type in AtomTypes:
            header = header + '  0.000\t' + '  0.000\t' + '  0.000\t' + '  0.000\t'
        header = header + '\n'

        with open(file_2, 'w') as f:
            f.write(header)
            for DataPoints in Time_and_MSD:
                f.write('{:7.3f}'.format(DataPoints[0])) # Time
                for ListPosition in range(0,len(AtomTypes)):
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2  ]   ))
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2+1][0]))
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2+1][1]))
                    f.write('\t{:7.3f}'.format(DataPoints[1 + ListPosition * 2+1][2]))
                f.write('\n')

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------!Conductivity!----------------------------------------------------------------# 

class cond:

    def conductivity(charge_carriers, volume, diff, temperature, hr):
        
        """
            Calculate the ionic conductivity.

        Args:
            charge_carriers (:py:attr:`float`): Number of charge carriers.
            volume (:py:attr:`float`): Average cell volume.
            diff (:py:attr:`float`): Charge Diffusion coefficient .
            temperature (:py:attr:`float`): Temperature.
            hr (:py:attr:`float`): Haven ratio if Tracer diffusion used otherwise put = 1
            Returns:
            conductivity (:py:attr:`float`): Ionic conductivity.
        """

        volume = volume * (10 ** -24)
        diff = diff * (10 ** -8)
        conc = charge_carriers / volume
        EV = ev ** 2
        constants = kb * temperature
        conductivity = ((diff * conc) * EV) / constants
        return conductivity * hr

class lammps_log():
    """
        Properties:- 

                    ['Density', 'PotEng', 'Pxx', 'Pxy', 'Pxz', 'Pyy', 'Pyz', 'Pzz', 'Step', 'Temp', 'TotEng', 'Volume'] etc.

                    It has capability to read all the fractional runs inside log file.

    """
    def tar_lammps(job,prop_1,prop_2):
        den = []
        step = []
        job_file_compressed = os.path.join(job.working_directory, [f for f in os.listdir(job.working_directory) if ".tar.bz2" in f][0])
        tar = tarfile.open(job_file_compressed, "r:bz2")
        for member in tar.getmembers():
            if 'log.lammps' == member.name:
                f = tar.extractfile(member)
                log = lammps_logfile.File(utf8reader(f))
                for i in range(log.get_num_partial_logs()):
                    x = log.get(prop_1, run_num=i)
                    y = log.get(prop_2, run_num=i)
                    den.append(y)
                    step.append(x)
        return step, den

    def tar_lammps_2(path,prop_1,prop_2,filename='log.lammps'):
        den = []
        step = []

        os.system('tar -xvf %s %s'%(path,filename))
        #os.system('mv log.lammps' + ' ' + path)
        #file = path + '/%s'%(filename)
        file = '%s'%(filename)

        log = lammps_logfile.File(file)
        for i in range(log.get_num_partial_logs()):
            x = log.get(prop_1, run_num=i)
            y = log.get(prop_2, run_num=i)
            den.append(y)
            step.append(x)
        return step, den

    def tar_lammps_3(job):
        job_file_compressed = os.path.join(job.working_directory, [f for f in os.listdir(job.working_directory) if ".tar.bz2" in f][0])
        tar = tarfile.open(job_file_compressed, "r:bz2")
        for member in tar.getmembers():
            if 'log.lammps' == member.name:
                f = tar.extractfile(member)
                log = lammps_logfile.File(utf8reader(f))
 
        return log

    def extract_file(job,filename='dump_nvt_prod.out'):
        
        job_file_compressed = os.path.join(job.working_directory, [f for f in os.listdir(job.working_directory) if ".tar.bz2" in f][0])
        os.system('tar -xvf %s %s'%(job_file_compressed,filename))
        os.system('mv dump_nvt_prod.out' + ' ' + job.working_directory)
        file = job.working_directory + '/%s'%(filename)

        return file

    def last_struct_ovito(file_2,path,filename='dump.out'):
        os.system('tar -xvf %s %s'%(path,filename))
        os.system('mv dump.out' + ' ' + path)
        file = path + '/%s'%(filename)
        pipeline = import_file(file)
        x = range(pipeline.source.num_frames)[-1]
        
        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)
        
        export_file(pipeline, file_2, Frame=x, 
            columns = ["Particle Type", "Position.X", "Position.Y", "Position.Z"])
        

    def extract_file_2(filename='dump_nvt_prod.out',path=True,path_line='',remove=False):

        if path==True:

            job_file_compressed = os.path.join(path_line, [f for f in os.listdir(path_line) if ".tar.bz2" in f][0])
            os.system('tar -xvf %s %s'%(job_file_compressed,filename))
            if remove == True:
                os.system('rm dump_nvt_prod.out')
            else:
                os.system('mv dump_nvt_prod.out' + ' ' + path_line)
            file = path_line + '/%s'%(filename)

        return file

        #tar = tarfile.open(job_file_compressed, "r:bz2")
        #for filename in tar.getmembers():
        #    shutil.unpack_archive(filename, job.working_directory)

    
    def pyiron_ovito_msd(job,filename='dump_nvt_prod.out',fil='cryt_msd_523k.txt'):#,fil_2='msd.txt'):

        file = lammps_log.extract_file(job=job,filename=filename)
        msd.tr_msd(file=file, file_2=fil)

    

class ovito_msd():

    def diff(job,filename='dump_nvt_prod.out',file_2='msd_na.txt',step=1.0, dump_write=100, step_ovito=10):

        job_file_compressed = os.path.join(job.working_directory, [f for f in os.listdir(job.working_directory) if ".tar.bz2" in f][0])
        os.system('tar -xvf %s %s'%(job_file_compressed,filename))
        os.system('mv dump_nvt_prod.out' + ' ' + job.working_directory)
        file = job.working_directory + '/%s'%(filename)

        #file = job.working_directory + '/%s'%(filename)

        # Data import 
        pipeline = import_file(file)

        # Simulation Parameters
        POTIM = step                  # Time step in fs
        NBLOCK = dump_write           # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
        # MSD Evaluation
        Stepsize = step_ovito

        pipeline.modifiers.append(SmoothTrajectoryModifier(minimum_image_convention=True,window_size=Stepsize))
        # Print the list of input particle type

        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)

        def modify1(frame: int, data: DataCollection):

            charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
            mass_dict = {"Zr": 91.224,
                            "Si": 28.086,
                            "P": 30.974,
                            "O": 15.999,
                            "Na": 22.99}
            charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
            mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
            ptypes = data.particles_.particle_types_
            for i in range(data.particles.count):
                charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
        pipeline.modifiers.append(modify1)

        # Displacement vectors

        pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

        # Calculate 'Relative Displacement' of host matrix

        def modify2(frame, data, type_id_Na = 2):
            
            #with respect to center of mass of host matrix                                                                                                                                                                                        
            ptypes = data.particles['Particle Type']
            displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
            displ = data.particles['Displacement']
            masses = data.particles['Mass'][(ptypes != type_id_Na)]
            
            #calculate center of mass displacement of host matrix                                                                                                                                                                                 
            total_mass_matrix = np.sum( masses)
            sum = 0.0
            for i in range(len(displ_matrix)):
                    sum +=displ_matrix[i] * masses[i]
            com_displ_host_matrix = sum/total_mass_matrix
            data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
            data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
    
        import functools
        pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

        def modify3(frame: int, data: DataCollection, type_id_Na = 2):

            ptype = data.particles.particle_type

            time = frame*NBLOCK*POTIM/1000

            bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
            bulk_displ = data.particles['Displacement'][(ptype == type_id_Na)]
        
            msd_x, msd_y, msd_z = np.mean(bulk_displ_rel**2, axis = 0)
            data.attributes["MSDx_bulk_com"] = msd_x  # com = center of mass /displacements are corrected for com shift of host matrix
            data.attributes["MSDy_bulk_com"] = msd_y  # MSD = tracer diffucivity
            data.attributes["MSDz_bulk_com"] = msd_z
            data.attributes["MSDav_bulk_com"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.mean(bulk_displ_rel, axis = 0)
            data.attributes["netMSDx_bulk_com"] = net_x**2  # netMSD = ionic diffusivity
            data.attributes["netMSDy_bulk_com"] = net_y**2
            data.attributes["netMSDz_bulk_com"] = net_z**2
            data.attributes["netMSDav_bulk_com"] = (net_x+net_y+net_z)**2

            msd_x, msd_y, msd_z = np.mean(bulk_displ**2, axis = 0)
            data.attributes["MSDx_bulk_raw"] = msd_x # raw = raw data/shift of host matrix is not considered
            data.attributes["MSDy_bulk_raw"] = msd_y 
            data.attributes["MSDz_bulk_raw"] = msd_z
            data.attributes["MSDav_bulk_raw"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.mean(bulk_displ, axis = 0)
            data.attributes["netMSDx_bulk_raw"] = net_x**2
            data.attributes["netMSDy_bulk_raw"] = net_y**2
            data.attributes["netMSDz_bulk_raw"] = net_z**2
            data.attributes["netMSDav_bulk_raw"] = (net_x+net_y+net_z)**2

            data.attributes["Timestep"] = time
   
        
        pipeline.modifiers.append(modify3)

        export_file(pipeline, file_2, "txt/attr",
        columns=["Timestep", "MSDx_bulk_raw", "MSDy_bulk_raw", "MSDz_bulk_raw","MSDav_bulk_raw","netMSDx_bulk_raw", "netMSDy_bulk_raw", "netMSDz_bulk_raw","netMSDav_bulk_raw","MSDx_bulk_com", "MSDy_bulk_com", "MSDz_bulk_com", "MSDav_bulk_com","netMSDx_bulk_com", "netMSDy_bulk_com", "netMSDz_bulk_com","netMSDav_bulk_com"],
        multiple_frames=True)
        
    def diff_ch(job,filename='dump_nvt_prod.out',file_2='msd_na.txt',step=1.0, dump_write=100, step_ovito=10):

        job_file_compressed = os.path.join(job.working_directory, [f for f in os.listdir(job.working_directory) if ".tar.bz2" in f][0])
        os.system('tar -xvf %s %s'%(job_file_compressed,filename))
        os.system('mv dump_nvt_prod.out' + ' ' + job.working_directory)
        file = job.working_directory + '/%s'%(filename)

        #file = job.working_directory + '/%s'%(filename)

        # Data import 
        pipeline = import_file(file)

        # Simulation Parameters
        POTIM = step                  # Time step in fs
        NBLOCK = dump_write           # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
        # MSD Evaluation
        Stepsize = step_ovito

        pipeline.modifiers.append(SmoothTrajectoryModifier(minimum_image_convention=True,window_size=Stepsize))
        # Print the list of input particle type

        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)

        def modify1(frame: int, data: DataCollection):

            charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
            mass_dict = {"Zr": 91.224,
                            "Si": 28.086,
                            "P": 30.974,
                            "O": 15.999,
                            "Na": 22.99}
            charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
            mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
            ptypes = data.particles_.particle_types_
            for i in range(data.particles.count):
                charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
        pipeline.modifiers.append(modify1)

        # Displacement vectors

        pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

        # Calculate 'Relative Displacement' of host matrix

        def modify2(frame, data, type_id_Na = 2):
            
            #with respect to center of mass of host matrix                                                                                                                                                                                        
            ptypes = data.particles['Particle Type']
            displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
            displ = data.particles['Displacement']
            masses = data.particles['Mass'][(ptypes != type_id_Na)]
            
            #calculate center of mass displacement of host matrix                                                                                                                                                                                 
            total_mass_matrix = np.sum( masses)
            sum = 0.0
            for i in range(len(displ_matrix)):
                    sum +=displ_matrix[i] * masses[i]
            com_displ_host_matrix = sum/total_mass_matrix
            data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
            data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
    
        import functools
        pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

        def modify3(frame: int, data: DataCollection, type_id_Na = 2):

            ptype = data.particles.particle_type

            time = frame*NBLOCK*POTIM/1000

            bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
            #bulk_displ = data.particles['Displacement'][(ptype == type_id_Na)]
        
            net_x, net_y, net_z = np.sum(bulk_displ_rel, axis = 0)
            data.attributes["netMSDx_bulk_com"] = np.mean(net_x**2) #/ len(bulk_displ_rel) #netMSD = ionic diffusivity
            data.attributes["netMSDy_bulk_com"] = np.mean(net_y**2) #/ len(bulk_displ_rel)
            data.attributes["netMSDz_bulk_com"] = np.mean(net_z**2) #/ len(bulk_displ_rel)
            data.attributes["netMSDav_bulk_com"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ_rel)
            
            data.attributes["Timestep"] = time
        
        pipeline.modifiers.append(modify3)

        export_file(pipeline, file_2, "txt/attr",
        columns=["Timestep","netMSDx_bulk_com", "netMSDy_bulk_com", "netMSDz_bulk_com","netMSDav_bulk_com"],
        multiple_frames=True)

    def diff_2(path_line='dump_nvt_prod.out',file_2='msd_na.txt',step=1.0, dump_write=100, step_ovito=10):

        file = path_line

        # Data import 
        pipeline = import_file(file)

        # Simulation Parameters
        POTIM = step                  # Time step in fs
        NBLOCK = dump_write            # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
        # MSD Evaluation
        Stepsize = step_ovito 

        # Print the list of input particle type
        pipeline.modifiers.append(SmoothTrajectoryModifier(minimum_image_convention=True,window_size=Stepsize))

        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)

        def modify1(frame: int, data: DataCollection):

            charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
            mass_dict = {"Zr": 91.224,
                            "Si": 28.086,
                            "P": 30.974,
                            "O": 15.999,
                            "Na": 22.99}
            charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
            mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
            ptypes = data.particles_.particle_types_
            for i in range(data.particles.count):
                charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
        pipeline.modifiers.append(modify1)
        
        # Cartesian to fractional coordinates 
        def cart_redu(frame, data):
            data = pipeline.compute()
            cartesian_positions     = data.particles.positions
            reduced_positions       = (data.cell.inverse[0:3,0:3] @ cartesian_positions.T).T + data.cell.inverse[0:3,3]
            #cartesian_positions_out = (cell[0:3,0:3] @ reduced_positions.T).T + cell[0:3,3]
            #assert np.allclose(cartesian_positions_out, cartesian_positions)
    
        pipeline.modifiers.append(cart_redu)

        # Displacement vectors

        pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

        # Calculate 'Relative Displacement' of host matrix

        def modify2(frame, data, type_id_Na = 2):
            
            #with respect to center of mass of host matrix                                                                                                                                                                                        
            ptypes = data.particles['Particle Type']
            displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
            displ = data.particles['Displacement']
            masses = data.particles['Mass']
            
            #calculate center of mass displacement of host matrix                                                                                                                                                                                 
            total_mass_matrix = np.sum( masses[ (ptypes != type_id_Na) ] )
            sum = 0.0
            for i in range(len(displ_matrix)):
                    sum +=displ_matrix[i] * masses[i]
            com_displ_host_matrix = sum/total_mass_matrix
            data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
            data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
    
        import functools
        pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

        def modify3(frame: int, data: DataCollection, type_id_Na = 2):

            ptype = data.particles.particle_type

            time = frame*NBLOCK*POTIM/1000

            bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
            bulk_displ = data.particles['Displacement'][(ptype == type_id_Na)]
        
            msd_x, msd_y, msd_z = np.mean(bulk_displ_rel**2, axis = 0)
            data.attributes["MSDx_bulk_com"] = msd_x  # com = center of mass /displacements are corrected for com shift of host matrix
            data.attributes["MSDy_bulk_com"] = msd_y  # MSD = tracer diffucivity
            data.attributes["MSDz_bulk_com"] = msd_z
            data.attributes["MSDav_bulk_com"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.sum(bulk_displ_rel, axis = 0)
            data.attributes["netMSDx_bulk_com"] = net_x**2 / len(bulk_displ_rel)  # netMSD = ionic diffusivity
            data.attributes["netMSDy_bulk_com"] = net_y**2 / len(bulk_displ_rel)
            data.attributes["netMSDz_bulk_com"] = net_z**2 / len(bulk_displ_rel)
            data.attributes["netMSDav_bulk_com"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ_rel)

            msd_x, msd_y, msd_z = np.mean(bulk_displ**2, axis = 0)
            data.attributes["MSDx_bulk_raw"] = msd_x # raw = raw data/shift of host matrix is not considered
            data.attributes["MSDy_bulk_raw"] = msd_y 
            data.attributes["MSDz_bulk_raw"] = msd_z
            data.attributes["MSDav_bulk_raw"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.sum(bulk_displ, axis = 0)
            data.attributes["netMSDx_bulk_raw"] = net_x**2 / len(bulk_displ)
            data.attributes["netMSDy_bulk_raw"] = net_y**2 / len(bulk_displ)
            data.attributes["netMSDz_bulk_raw"] = net_z**2 / len(bulk_displ)
            data.attributes["netMSDav_bulk_raw"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ)

            data.attributes["Timestep"] = time
   
        
        pipeline.modifiers.append(modify3)

        export_file(pipeline, file_2, "txt/attr",
        columns=["Timestep", "MSDx_bulk_raw", "MSDy_bulk_raw", "MSDz_bulk_raw","MSDav_bulk_raw","netMSDx_bulk_raw", "netMSDy_bulk_raw", "netMSDz_bulk_raw","netMSDav_bulk_raw","MSDx_bulk_com", "MSDy_bulk_com", "MSDz_bulk_com", "MSDav_bulk_com","netMSDx_bulk_com", "netMSDy_bulk_com", "netMSDz_bulk_com","netMSDav_bulk_com"],
        multiple_frames=True)

    def diff_3(path_line='dump_nvt_prod.out',file_2='msd_na.txt',step=1.0, dump_write=100, step_ovito=10):

        file = path_line

        # Data import 
        pipeline = import_file(file)

        # Simulation Parameters
        POTIM = step                  # Time step in fs
        NBLOCK = dump_write            # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
        # MSD Evaluation
        Stepsize = step_ovito 

        # Print the list of input particle type
        pipeline.modifiers.append(SmoothTrajectoryModifier(minimum_image_convention=True,window_size=Stepsize))

        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)

        def modify1(frame: int, data: DataCollection):

            charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
            mass_dict = {"Zr": 91.224,
                            "Si": 28.086,
                            "P": 30.974,
                            "O": 15.999,
                            "Na": 22.99}
            charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
            mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
            ptypes = data.particles_.particle_types_
            for i in range(data.particles.count):
                charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
        pipeline.modifiers.append(modify1)

        # Rotation of the box to 30 degree clockwise to extraxt Dyy 
        # Affine transformation:
        pipeline.modifiers.append(AffineTransformationModifier(transformation = [[0.8660254037844387, 0.49999999999999994, 0.0, 0.0], [-0.49999999999999994, 0.8660254037844387, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]))

        # Displacement vectors

        pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

        # Calculate 'Relative Displacement' of host matrix

        def modify2(frame, data, type_id_Na = 2):
            
            #with respect to center of mass of host matrix                                                                                                                                                                                        
            ptypes = data.particles['Particle Type']
            displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
            displ = data.particles['Displacement']
            masses = data.particles['Mass']
            
            #calculate center of mass displacement of host matrix                                                                                                                                                                                 
            total_mass_matrix = np.sum( masses[ (ptypes != type_id_Na) ] )
            sum = 0.0
            for i in range(len(displ_matrix)):
                    sum +=displ_matrix[i] * masses[i]
            com_displ_host_matrix = sum/total_mass_matrix
            data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
            data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
    
        import functools
        pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

        def modify3(frame: int, data: DataCollection, type_id_Na = 2):

            ptype = data.particles.particle_type

            time = frame*NBLOCK*POTIM/1000

            bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
            bulk_displ = data.particles['Displacement'][(ptype == type_id_Na)]
        
            msd_x, msd_y, msd_z = np.mean(bulk_displ_rel**2, axis = 0)
            data.attributes["MSDx_bulk_com"] = msd_x  # com = center of mass /displacements are corrected for com shift of host matrix
            data.attributes["MSDy_bulk_com"] = msd_y  # MSD = tracer diffucivity
            data.attributes["MSDz_bulk_com"] = msd_z
            data.attributes["MSDav_bulk_com"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.sum(bulk_displ_rel, axis = 0)
            data.attributes["netMSDx_bulk_com"] = net_x**2 / len(bulk_displ_rel) # netMSD = ionic diffusivity
            data.attributes["netMSDy_bulk_com"] = net_y**2 / len(bulk_displ_rel)
            data.attributes["netMSDz_bulk_com"] = net_z**2 / len(bulk_displ_rel)
            data.attributes["netMSDav_bulk_com"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ_rel)

            msd_x, msd_y, msd_z = np.mean(bulk_displ**2, axis = 0)
            data.attributes["MSDx_bulk_raw"] = msd_x # raw = raw data/shift of host matrix is not considered
            data.attributes["MSDy_bulk_raw"] = msd_y 
            data.attributes["MSDz_bulk_raw"] = msd_z
            data.attributes["MSDav_bulk_raw"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.sum(bulk_displ, axis = 0)
            data.attributes["netMSDx_bulk_raw"] = net_x**2 / len(bulk_displ)
            data.attributes["netMSDy_bulk_raw"] = net_y**2 / len(bulk_displ)
            data.attributes["netMSDz_bulk_raw"] = net_z**2 / len(bulk_displ)
            data.attributes["netMSDav_bulk_raw"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ)

            data.attributes["Timestep"] = time
   
        
        pipeline.modifiers.append(modify3)

        export_file(pipeline, file_2, "txt/attr",
        columns=["Timestep", "MSDx_bulk_raw", "MSDy_bulk_raw", "MSDz_bulk_raw","MSDav_bulk_raw","netMSDx_bulk_raw", "netMSDy_bulk_raw", "netMSDz_bulk_raw","netMSDav_bulk_raw","MSDx_bulk_com", "MSDy_bulk_com", "MSDz_bulk_com", "MSDav_bulk_com","netMSDx_bulk_com", "netMSDy_bulk_com", "netMSDz_bulk_com","netMSDav_bulk_com"],
        multiple_frames=True)

    def diff_4(path_line='dump_nvt_prod.out',file_2='msd_na.txt',step=1.0, dump_write=100, step_ovito=10):

        file = path_line

        # Data import 
        pipeline = import_file(file)

        # Simulation Parameters
        POTIM = step                  # Time step in fs
        NBLOCK = dump_write            # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
        # MSD Evaluation
        Stepsize = step_ovito 

        # Print the list of input particle type
        pipeline.modifiers.append(SmoothTrajectoryModifier(minimum_image_convention=True,window_size=Stepsize))
        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)

        def modify1(frame: int, data: DataCollection):

            charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
            mass_dict = {"Zr": 91.224,
                            "Si": 28.086,
                            "P": 30.974,
                            "O": 15.999,
                            "Na": 22.99}
            charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
            mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
            ptypes = data.particles_.particle_types_
            for i in range(data.particles.count):
                charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
        pipeline.modifiers.append(modify1)
        
        # Cartesian to fractional coordinates

        or_cell = []
        or_cell_2 = []
    
        def trans(frame, data):
            tm = data.cell.inverse[...] # Matrix which will tranform from cartesian to fractional, indirect
            tm_n = data.cell[...]
            #print(tm_n[0:3,0:3])
            or_cell.extend(tm_n)
            or_cell_2.extend(data.cell[0:3,3])
            # Execute AffineTransformationModifier as a sub-operation:
            #print((or_cell_2))
            data.apply(AffineTransformationModifier(transformation = tm))

        pipeline.modifiers.append(trans)


        pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

        # Calculate 'Relative Displacement' of host matrix

        def modify2(frame, data, type_id_Na = 2):
            
            #with respect to center of mass of host matrix                                                                                                                                                                                      
            ptypes = data.particles['Particle Type']
            #displ_matrix = reduced_positions['Displacement'][ (ptypes != type_id_Na) ]
            displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
            displ = data.particles['Displacement']
            masses = data.particles['Mass']                                                                                                                                                                           
            total_mass_matrix = np.sum( masses[ (ptypes != type_id_Na) ] )
            sum = 0.0
            for i in range(len(displ_matrix)):
                    sum +=displ_matrix[i] * masses[i]
            com_displ_host_matrix = sum/total_mass_matrix
            data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
            data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
    
        import functools
        pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

        from numpy.linalg import inv
        
        def modify3(frame: int, data: DataCollection, type_id_Na = 2):

            ptype = data.particles.particle_type

            time = frame*NBLOCK*POTIM/1000
            A = np.vstack([or_cell[0], or_cell[1], or_cell[2]])
            B = [or_cell_2]
            C = B[0][len(B[0])-3:len(B[0])]

            # Convert from fractional coordinate diplacement to cartesian 
            bulk_displ_rel = (A[0:3,0:3] @ data.particles['Relative Displacement'][(ptype == type_id_Na)].T).T + C
            bulk_displ = (A[0:3,0:3] @ data.particles['Displacement'][(ptype == type_id_Na)].T).T + C
        

            msd_x, msd_y, msd_z = np.mean(bulk_displ_rel**2, axis = 0)
            data.attributes["MSDx_bulk_com"] = msd_x  # com = center of mass /displacements are corrected for com shift of host matrix
            data.attributes["MSDy_bulk_com"] = msd_y  # MSD = tracer diffucivity
            data.attributes["MSDz_bulk_com"] = msd_z
            data.attributes["MSDav_bulk_com"] = msd_x+msd_y+msd_z 

            net_x, net_y, net_z = np.sum(bulk_displ_rel, axis = 0)
            data.attributes["netMSDx_bulk_com"] = net_x**2 / len(bulk_displ_rel)  # netMSD = ionic diffusivity
            data.attributes["netMSDy_bulk_com"] = net_y**2 / len(bulk_displ_rel)
            data.attributes["netMSDz_bulk_com"] = net_z**2 / len(bulk_displ_rel)
            data.attributes["netMSDav_bulk_com"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ_rel)

            msd_x, msd_y, msd_z = np.mean(bulk_displ**2, axis = 0)
            data.attributes["MSDx_bulk_raw"] = msd_x # raw = raw data/shift of host matrix is not considered
            data.attributes["MSDy_bulk_raw"] = msd_y 
            data.attributes["MSDz_bulk_raw"] = msd_z
            data.attributes["MSDav_bulk_raw"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.sum(bulk_displ, axis = 0)
            data.attributes["netMSDx_bulk_raw"] = net_x**2 / len(bulk_displ)
            data.attributes["netMSDy_bulk_raw"] = net_y**2 / len(bulk_displ)
            data.attributes["netMSDz_bulk_raw"] = net_z**2 / len(bulk_displ)
            data.attributes["netMSDav_bulk_raw"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ)

            data.attributes["Timestep"] = time
        
        pipeline.modifiers.append(modify3)

        export_file(pipeline, file_2, "txt/attr",
        columns=["Timestep", "MSDx_bulk_raw", "MSDy_bulk_raw", "MSDz_bulk_raw","MSDav_bulk_raw","netMSDx_bulk_raw", "netMSDy_bulk_raw", "netMSDz_bulk_raw","netMSDav_bulk_raw","MSDx_bulk_com", "MSDy_bulk_com", "MSDz_bulk_com", "MSDav_bulk_com","netMSDx_bulk_com", "netMSDy_bulk_com", "netMSDz_bulk_com","netMSDav_bulk_com"],
        multiple_frames=True)


    def regional_diff_twist(job,filename='dump_nvt_prod.out',file_2='msd_na.txt',w=10):

        file = job.working_directory + '/%s'%(filename)

        # Data import 
        pipeline = import_file(file)

        # Simulation Parameters
        POTIM = step = 1.00                 # Time step in fs
        NBLOCK = dump_write = 100           # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
        # MSD Evaluation
        Stepsize = step_ovito = 10

        # Print the list of input particle type
        pipeline.modifiers.append(SmoothTrajectoryModifier(minimum_image_convention=True,window_size=Stepsize))

        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)

        def modify1(frame: int, data: DataCollection):

            charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
            mass_dict = {"Zr": 91.224,
                            "Si": 28.086,
                            "P": 30.974,
                            "O": 15.999,
                            "Na": 22.99}
            charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
            mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
            ptypes = data.particles_.particle_types_
            for i in range(data.particles.count):
                charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
        pipeline.modifiers.append(modify1)

        # Displacement vectors

        pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

        # Calculate 'Relative Displacement' of host matrix

        def modify2(frame, data, type_id_Na = 2):
            
            #with respect to center of mass of host matrix                                                                                                                                                                                        
            ptypes = data.particles['Particle Type']
            displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
            displ = data.particles['Displacement']
            masses = data.particles['Mass']
            
            #calculate center of mass displacement of host matrix                                                                                                                                                                                 
            total_mass_matrix = np.sum( masses[ (ptypes != type_id_Na) ] )
            sum = 0.0
            for i in range(len(displ_matrix)):
                    sum +=displ_matrix[i] * masses[i]
            com_displ_host_matrix = sum/total_mass_matrix
            data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
            data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
    
        import functools
        pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

        def modify3(frame: int, data: DataCollection, type_id_Na = 2):

            #Select GB: Choice has to be double-checked in Ovito!
            
            pos = data.particles.positions
            Selection = data.particles_.create_property("Selection")
            Selection[(pos[:,2] < data.cell[2][3] + w) ] = 1
            Selection[(pos[:,2] > data.cell[2][3] + data.cell[2][2]/2. - w) & (pos[:,2] < data.cell[2][3] + data.cell[2][2]/2. + w) ] = 1
            Selection[(pos[:,2] > data.cell[2][3] + data.cell[2][2] - w) ] = 1
        
            data.apply(FreezePropertyModifier(source_property = 'Selection', 
                                      destination_property = 'Select0',
                                      freeze_at = 0))

            ptype = data.particles.particle_type
            #only Na atoms
            Selection = data.particles["Select0"]
            #print(type(Selection))
            time = frame*NBLOCK*POTIM/1000
            #print(time)

            gb_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na) & (Selection)]
            grain_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na) & (Selection != 1)]
            gb_displ = data.particles['Displacement'][(ptype == type_id_Na) & (Selection)]
            grain_displ = data.particles['Displacement'][(ptype == type_id_Na) & (Selection != 1)]
        
            msd_x, msd_y, msd_z = np.mean(gb_displ_rel**2, axis = 0)
            data.attributes["MSDx_gb_com"] = msd_x  # com = center of mass /displacements are corrected for com shift of host matrix
            data.attributes["MSDy_gb_com"] = msd_y  # MSD = tracer diffucivity
            data.attributes["MSDz_gb_com"] = msd_z
            data.attributes["MSDav_gb_com"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.mean(gb_displ_rel, axis = 0)
            data.attributes["netMSDx_gb_com"] = net_x**2  # netMSD = ionic diffusivity
            data.attributes["netMSDy_gb_com"] = net_y**2
            data.attributes["netMSDz_gb_com"] = net_z**2
            data.attributes["netMSDav_gb_com"] = (net_x+net_y+net_z)**2

            msd_x, msd_y, msd_z = np.mean(gb_displ**2, axis = 0)
            data.attributes["MSDx_gb_raw"] = msd_x # raw = raw data/shift of host matrix is not considered
            data.attributes["MSDy_gb_raw"] = msd_y 
            data.attributes["MSDz_gb_raw"] = msd_z
            data.attributes["MSDav_gb_raw"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.mean(gb_displ, axis = 0)
            data.attributes["netMSDx_gb_raw"] = net_x**2
            data.attributes["netMSDy_gb_raw"] = net_y**2
            data.attributes["netMSDz_gb_raw"] = net_z**2
            data.attributes["netMSDav_gb_raw"] = (net_x+net_y+net_z)**2

        
            msd_x, msd_y, msd_z = np.mean(grain_displ_rel**2, axis = 0)
            data.attributes["MSDx_grain_com"] = msd_x  #com = center of mass /displacements are corrected for com shift of host matrix
            data.attributes["MSDy_grain_com"] = msd_y
            data.attributes["MSDz_grain_com"] = msd_z
            data.attributes["MSDav_grain_com"] = msd_x+msd_y+msd_z
            

            net_x, net_y, net_z = np.mean(grain_displ_rel, axis = 0)
            data.attributes["netMSDx_grain_com"] = net_x**2
            data.attributes["netMSDy_grain_com"] = net_y**2
            data.attributes["netMSDz_grain_com"] = net_z**2
            data.attributes["netMSDav_grain_com"] = (net_x+net_y+net_z)**2

            msd_x, msd_y, msd_z = np.mean(grain_displ**2, axis = 0)
            data.attributes["MSDx_grain_raw"] = msd_x # raw = raw data/shift of host matrix is not considered
            data.attributes["MSDy_grain_raw"] = msd_y 
            data.attributes["MSDz_grain_raw"] = msd_z
            data.attributes["MSDav_grain_raw"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.mean(grain_displ, axis = 0)
            data.attributes["netMSDx_grain_raw"] = net_x**2
            data.attributes["netMSDy_grain_raw"] = net_y**2
            data.attributes["netMSDz_grain_raw"] = net_z**2
            data.attributes["netMSDav_grain_raw"] = (net_x+net_y+net_z)**2

            data.attributes["Timestep"] = time
   
        
        pipeline.modifiers.append(modify3)

        export_file(pipeline, file_2, "txt/attr",
        columns=["Timestep", "MSDx_gb_raw", "MSDy_gb_raw", "MSDz_gb_raw", "MSDav_gb_raw", "netMSDx_gb_raw", "netMSDy_gb_raw", "netMSDz_gb_raw","netMSDav_gb_raw", "MSDx_gb_com", "MSDy_gb_com", "MSDz_gb_com","MSDav_gb_com", "netMSDx_gb_com", "netMSDy_gb_com", "netMSDz_gb_com","netMSDav_gb_com",
                "MSDx_grain_raw", "MSDy_grain_raw", "MSDz_grain_raw","MSDav_grain_raw", "netMSDx_grain_raw", "netMSDy_grain_raw", "netMSDz_grain_raw","netMSDav_grain_raw", "MSDx_grain_com", "MSDy_grain_com", "MSDz_grain_com","MSDav_grain_com", "netMSDx_grain_com", "netMSDy_grain_com", "netMSDz_grain_com","netMSDav_grain_com"],
        multiple_frames=True)


    def regional_diff_tilt(path_line='dump_nvt_prod.out',file_2='msd_na.txt',w=10):

        file = path_line

        # Data import 
        pipeline = import_file(file)

        # Simulation Parameters
        POTIM = step = 1.00                 # Time step in fs
        NBLOCK = dump_write = 100           # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
        # MSD Evaluation
        Stepsize = step_ovito = 10

        # Trajectory smoothning modifier 
        pipeline.modifiers.append(SmoothTrajectoryModifier(minimum_image_convention=True,window_size=Stepsize))

        # Print the list of input particle type

        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)

        def modify1(frame: int, data: DataCollection):

            charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
            mass_dict = {"Zr": 91.224,
                            "Si": 28.086,
                            "P": 30.974,
                            "O": 15.999,
                            "Na": 22.99}
            charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
            mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
            ptypes = data.particles_.particle_types_
            for i in range(data.particles.count):
                charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
        pipeline.modifiers.append(modify1)

        # Displacement vectors

        pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

        # Calculate 'Relative Displacement' of host matrix

        def modify2(frame, data, type_id_Na = 2):
            
            #with respect to center of mass of host matrix                                                                                                                                                                                        
            ptypes = data.particles['Particle Type']
            displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
            displ = data.particles['Displacement']
            masses = data.particles['Mass']
            
            #calculate center of mass displacement of host matrix                                                                                                                                                                                 
            total_mass_matrix = np.sum( masses[ (ptypes != type_id_Na) ] )
            sum = 0.0
            for i in range(len(displ_matrix)):
                    sum +=displ_matrix[i] * masses[i]
            com_displ_host_matrix = sum/total_mass_matrix
            data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
            data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
    
        import functools
        pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

        def modify3(frame: int, data: DataCollection, type_id_Na = 2):

            #Select GB: Choice has to be double-checked in Ovito!
            
            pos = data.particles.positions
            Selection = data.particles_.create_property("Selection")
            Selection[(pos[:,0] < data.cell[0][3] + w) ] = 1
            Selection[(pos[:,0] > data.cell[0][3] + data.cell[0][0]/2. - w) & (pos[:,0] < data.cell[0][3] + data.cell[0][0]/2. + w) ] = 1
            Selection[(pos[:,0] > data.cell[0][3] + data.cell[0][0] - w) ] = 1
        
            data.apply(FreezePropertyModifier(source_property = 'Selection', 
                                      destination_property = 'Select0',
                                      freeze_at = 0))

            ptype = data.particles.particle_type
            #only Na atoms
            Selection = data.particles["Select0"]
            #print(type(Selection))
            time = frame*NBLOCK*POTIM/1000
            #print(time)

            gb_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na) & (Selection)]
            grain_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na) & (Selection != 1)]
            gb_displ = data.particles['Displacement'][(ptype == type_id_Na) & (Selection)]
            grain_displ = data.particles['Displacement'][(ptype == type_id_Na) & (Selection != 1)]
        
            msd_x, msd_y, msd_z = np.mean(gb_displ_rel**2, axis = 0)
            data.attributes["MSDx_gb_com"] = msd_x  # com = center of mass /displacements are corrected for com shift of host matrix
            data.attributes["MSDy_gb_com"] = msd_y  # MSD = tracer diffucivity
            data.attributes["MSDz_gb_com"] = msd_z
            data.attributes["MSDav_gb_com"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.mean(gb_displ_rel, axis = 0)
            data.attributes["netMSDx_gb_com"] = net_x**2  # netMSD = ionic diffusivity
            data.attributes["netMSDy_gb_com"] = net_y**2
            data.attributes["netMSDz_gb_com"] = net_z**2
            data.attributes["netMSDav_gb_com"] = (net_x+net_y+net_z)**2

            msd_x, msd_y, msd_z = np.mean(gb_displ**2, axis = 0)
            data.attributes["MSDx_gb_raw"] = msd_x # raw = raw data/shift of host matrix is not considered
            data.attributes["MSDy_gb_raw"] = msd_y 
            data.attributes["MSDz_gb_raw"] = msd_z
            data.attributes["MSDav_gb_raw"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.mean(gb_displ, axis = 0)
            data.attributes["netMSDx_gb_raw"] = net_x**2
            data.attributes["netMSDy_gb_raw"] = net_y**2
            data.attributes["netMSDz_gb_raw"] = net_z**2
            data.attributes["netMSDav_gb_raw"] = (net_x+net_y+net_z)**2

        
            msd_x, msd_y, msd_z = np.mean(grain_displ_rel**2, axis = 0)
            data.attributes["MSDx_grain_com"] = msd_x  #com = center of mass /displacements are corrected for com shift of host matrix
            data.attributes["MSDy_grain_com"] = msd_y
            data.attributes["MSDz_grain_com"] = msd_z
            data.attributes["MSDav_grain_com"] = msd_x+msd_y+msd_z
            

            net_x, net_y, net_z = np.mean(grain_displ_rel, axis = 0)
            data.attributes["netMSDx_grain_com"] = net_x**2
            data.attributes["netMSDy_grain_com"] = net_y**2
            data.attributes["netMSDz_grain_com"] = net_z**2
            data.attributes["netMSDav_grain_com"] = (net_x+net_y+net_z)**2

            msd_x, msd_y, msd_z = np.mean(grain_displ**2, axis = 0)
            data.attributes["MSDx_grain_raw"] = msd_x # raw = raw data/shift of host matrix is not considered
            data.attributes["MSDy_grain_raw"] = msd_y 
            data.attributes["MSDz_grain_raw"] = msd_z
            data.attributes["MSDav_grain_raw"] = msd_x+msd_y+msd_z

            net_x, net_y, net_z = np.mean(grain_displ, axis = 0)
            data.attributes["netMSDx_grain_raw"] = net_x**2
            data.attributes["netMSDy_grain_raw"] = net_y**2
            data.attributes["netMSDz_grain_raw"] = net_z**2
            data.attributes["netMSDav_grain_raw"] = (net_x+net_y+net_z)**2

            data.attributes["Timestep"] = time
   
        
        pipeline.modifiers.append(modify3)

        export_file(pipeline, file_2, "txt/attr",
        columns=["Timestep", "MSDx_gb_raw", "MSDy_gb_raw", "MSDz_gb_raw", "MSDav_gb_raw", "netMSDx_gb_raw", "netMSDy_gb_raw", "netMSDz_gb_raw","netMSDav_gb_raw", "MSDx_gb_com", "MSDy_gb_com", "MSDz_gb_com","MSDav_gb_com", "netMSDx_gb_com", "netMSDy_gb_com", "netMSDz_gb_com","netMSDav_gb_com",
                "MSDx_grain_raw", "MSDy_grain_raw", "MSDz_grain_raw","MSDav_grain_raw", "netMSDx_grain_raw", "netMSDy_grain_raw", "netMSDz_grain_raw","netMSDav_grain_raw", "MSDx_grain_com", "MSDy_grain_com", "MSDz_grain_com","MSDav_grain_com", "netMSDx_grain_com", "netMSDy_grain_com", "netMSDz_grain_com","netMSDav_grain_com"],
        multiple_frames=True)
        
class new_ovito:
    
    def tracer_diff(path_line='dump_nvt_prod.out',file_2='msd_na.txt',step=1.0, dump_write=100, step_ovito=10):

        file = path_line

        # Data import 
        pipeline = import_file(file)

        # Simulation Parameters
        POTIM = step                  # Time step in fs
        NBLOCK = dump_write            # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
        # MSD Evaluation
        Stepsize = step_ovito 

        # Print the list of input particle type
        pipeline.modifiers.append(SmoothTrajectoryModifier(minimum_image_convention=True,window_size=Stepsize))

        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)

        def modify1(frame: int, data: DataCollection):

            charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
            mass_dict = {"Zr": 91.224,
                            "Si": 28.086,
                            "P": 30.974,
                            "O": 15.999,
                            "Na": 22.99}
            charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
            mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
            ptypes = data.particles_.particle_types_
            for i in range(data.particles.count):
                charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
        pipeline.modifiers.append(modify1)

        # Displacement vectors

        pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

        # Calculate 'Relative Displacement' of host matrix

        def modify2(frame, data, type_id_Na = 2):
            
            #with respect to center of mass of host matrix                                                                                                                                                                                        
            ptypes = data.particles['Particle Type']
            displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
            displ = data.particles['Displacement']
            masses = data.particles['Mass'][(ptypes != type_id_Na)]
            
            #calculate center of mass displacement of host matrix                                                                                                                                                                                 
            total_mass_matrix = np.sum( masses)
            sum = 0.0
            for i in range(len(displ_matrix)):
                    sum +=displ_matrix[i] * masses[i]
            com_displ_host_matrix = sum/total_mass_matrix
            data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
            data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
    
        import functools
        pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

        def modify3(frame: int, data: DataCollection, type_id_Na = 2):

            ptype = data.particles.particle_type

            time = frame*NBLOCK*POTIM/1000

            pos = data.particles['Position'][(ptype == type_id_Na)]
            #print(len(pos[:,0]))
            bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
            bulk_displ = data.particles['Displacement'][(ptype == type_id_Na)]
            #print(bulk_displ_rel)
            # Note:- MSD are the absolute displacement suqared, therefore we need absolute displacement 
            #        The reason to have negative MSD in xy plane  
            msd_x, msd_y, msd_z = np.mean(np.abs(bulk_displ_rel)**2, axis = 0)
            data.attributes["MSDx_bulk_com"] = msd_x  # com = center of mass /displacements are corrected for com shift of host matrix
            data.attributes["MSDy_bulk_com"] = msd_y  # MSD = tracer diffucivity
            data.attributes["MSDz_bulk_com"] = msd_z
            data.attributes["MSDav_bulk_com"] = msd_x+msd_y+msd_z
        
            # The follwoing code implement cros_displacement algo to get offdiagonal compoenets of MSD 
            # Basically I am caulcating the varienece 
            #sum_of_pos_xy = np.mean((np.abs(bulk_displ_rel[:,0]) + np.abs(bulk_displ_rel[:,1]))**2)
            #sum_of_pos_xz = np.mean((np.abs(bulk_displ_rel[:,0]) + np.abs(bulk_displ_rel[:,2]))**2)
            #sum_of_pos_yz = np.mean((np.abs(bulk_displ_rel[:,1]) + np.abs(bulk_displ_rel[:,2]))**2)
        
            #cd_xy = cd[0][1]
            #cd_xz = cd[0][2]
            #cd_yz = cd[1][2]
        
            #print(cd_xy)
        
            #data.attributes["MSDxy_bulk_com"] = 0.5*(sum_of_pos_xy - msd_x - msd_y) #np.mean(cd_xy**2) #np.mean(np.array(cd_xy)**2, axis=1) 
            #data.attributes["MSDxz_bulk_com"] = 0.5*(sum_of_pos_xz - msd_x - msd_z) #np.mean(cd_xz**2) #np.mean(np.array(cd_xz)**2, axis=1) 
            #data.attributes["MSDyz_bulk_com"] = 0.5*(sum_of_pos_yz - msd_y - msd_z) #np.mean(cd_yz**2) #np.mean(np.array(cd_yz)**2, axis=1)

            #net_x, net_y, net_z = np.sum(np.abs(bulk_displ_rel), axis = 0)
            #data.attributes["netMSDx_bulk_com"] = net_x**2 / len(bulk_displ_rel) #netMSD = ionic diffusivity
            #data.attributes["netMSDy_bulk_com"] = net_y**2 / len(bulk_displ_rel)
            #data.attributes["netMSDz_bulk_com"] = net_z**2 / len(bulk_displ_rel)
            #data.attributes["netMSDav_bulk_com"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ_rel)
        
            # The follwoing code implement cros_displacement algo to get offdiagonal compoenets of MSD 
            #sum_of_pos_xy_n = ((np.sum(np.abs(bulk_displ_rel[:,0])) + np.abs(np.sum(bulk_displ_rel[:,1])))**2) / len(bulk_displ_rel)
            #sum_of_pos_xz_n = ((np.sum(np.abs(bulk_displ_rel[:,0])) + np.abs(np.sum(bulk_displ_rel[:,2])))**2) / len(bulk_displ_rel)
            #sum_of_pos_yz_n = ((np.sum(np.abs(bulk_displ_rel[:,1])) + np.abs(np.sum(bulk_displ_rel[:,2])))**2) / len(bulk_displ_rel)
        
            #data.attributes["netMSDxy_bulk_com"] = 0.5*(sum_of_pos_xy_n - net_x**2/len(np.abs(bulk_displ_rel)) - net_y**2/len(np.abs(bulk_displ_rel))) #np.mean(cd_xy**2) #np.mean(np.array(cd_xy)**2, axis=1) 
            #data.attributes["netMSDxz_bulk_com"] = 0.5*(sum_of_pos_xz_n - net_x**2/len(np.abs(bulk_displ_rel)) - net_z**2/len(np.abs(bulk_displ_rel))) #np.mean(cd_xz**2) #np.mean(np.array(cd_xz)**2, axis=1) 
            #data.attributes["netMSDyz_bulk_com"] = 0.5*(sum_of_pos_yz_n - net_y**2/len(np.abs(bulk_displ_rel)) - net_z**2/len(np.abs(bulk_displ_rel)))
        
            data.attributes["Timestep"] = time
     
        pipeline.modifiers.append(modify3)

        export_file(pipeline, file_2, "txt/attr",
        columns=["Timestep", "MSDx_bulk_com", "MSDy_bulk_com", "MSDz_bulk_com", "MSDav_bulk_com"],
        multiple_frames=True)    
    
    def charge_diff(path_line='dump_nvt_prod.out',file_2='msd_na.txt',step=1.0, dump_write=100, step_ovito=10):

        file = path_line

        # Data import 
        pipeline = import_file(file)

        # Simulation Parameters
        POTIM = step                  # Time step in fs
        NBLOCK = dump_write            # Each NBLOCK step is written to XDATCAR; for lammps -> after what amount of steps do you write out a dump file?
        
        # MSD Evaluation
        Stepsize = step_ovito 

        # Print the list of input particle type
        pipeline.modifiers.append(SmoothTrajectoryModifier(minimum_image_convention=True,window_size=Stepsize))

        # Setting up the atom type 
        pipeline.add_to_scene()
        def setup_particle_types(frame, data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "O"
            types.type_by_id_(2).name = "Na"
            types.type_by_id_(3).name = "Zr"
            types.type_by_id_(4).name = "Si"
            types.type_by_id_(5).name = "P"
        pipeline.modifiers.append(setup_particle_types)

        def modify1(frame: int, data: DataCollection):

            charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
            mass_dict = {"Zr": 91.224,
                            "Si": 28.086,
                            "P": 30.974,
                            "O": 15.999,
                            "Na": 22.99}
            charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
            mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
            ptypes = data.particles_.particle_types_
            for i in range(data.particles.count):
                charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
        pipeline.modifiers.append(modify1)

        # Displacement vectors

        pipeline.modifiers.append(CalculateDisplacementsModifier(minimum_image_convention = False))

        # Calculate 'Relative Displacement' of host matrix

        def modify2(frame, data, type_id_Na = 2):
            
            #with respect to center of mass of host matrix                                                                                                                                                                                        
            ptypes = data.particles['Particle Type']
            displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
            displ = data.particles['Displacement']
            masses = data.particles['Mass'][(ptypes != type_id_Na)]
            
            #calculate center of mass displacement of host matrix                                                                                                                                                                                 
            total_mass_matrix = np.sum( masses)
            sum = 0.0
            for i in range(len(displ_matrix)):
                    sum +=displ_matrix[i] * masses[i]
            com_displ_host_matrix = sum/total_mass_matrix
            data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
            data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
    
        import functools
        pipeline.modifiers.append(functools.partial(modify2, type_id_Na = 2))

        def modify3(frame: int, data: DataCollection, type_id_Na = 2):

            ptype = data.particles.particle_type

            time = frame*NBLOCK*POTIM/1000

            pos = data.particles['Position'][(ptype == type_id_Na)]
            #print(len(pos[:,0]))
            bulk_displ_rel = data.particles['Relative Displacement'][(ptype == type_id_Na)]
            bulk_displ = data.particles['Displacement'][(ptype == type_id_Na)]
            #print(bulk_displ_rel)
            # Note:- MSD are the absolute displacement suqared, therefore we need absolute displacement 
            #        The reason to have negative MSD in xy plane  
            #msd_x, msd_y, msd_z = np.mean(np.abs(bulk_displ_rel)**2, axis = 0)
            #data.attributes["MSDx_bulk_com"] = msd_x  # com = center of mass /displacements are corrected for com shift of host matrix
            #data.attributes["MSDy_bulk_com"] = msd_y  # MSD = tracer diffucivity
            #data.attributes["MSDz_bulk_com"] = msd_z
            #data.attributes["MSDav_bulk_com"] = msd_x+msd_y+msd_z
        
            # The follwoing code implement cros_displacement algo to get offdiagonal compoenets of MSD 
            # Basically I am caulcating the varienece 
            #sum_of_pos_xy = np.mean((np.abs(bulk_displ_rel[:,0]) + np.abs(bulk_displ_rel[:,1]))**2)
            #sum_of_pos_xz = np.mean((np.abs(bulk_displ_rel[:,0]) + np.abs(bulk_displ_rel[:,2]))**2)
            #sum_of_pos_yz = np.mean((np.abs(bulk_displ_rel[:,1]) + np.abs(bulk_displ_rel[:,2]))**2)
        
            #cd_xy = cd[0][1]
            #cd_xz = cd[0][2]
            #cd_yz = cd[1][2]
        
            #print(cd_xy)
        
            #data.attributes["MSDxy_bulk_com"] = 0.5*(sum_of_pos_xy - msd_x - msd_y) #np.mean(cd_xy**2) #np.mean(np.array(cd_xy)**2, axis=1) 
            #data.attributes["MSDxz_bulk_com"] = 0.5*(sum_of_pos_xz - msd_x - msd_z) #np.mean(cd_xz**2) #np.mean(np.array(cd_xz)**2, axis=1) 
            #data.attributes["MSDyz_bulk_com"] = 0.5*(sum_of_pos_yz - msd_y - msd_z) #np.mean(cd_yz**2) #np.mean(np.array(cd_yz)**2, axis=1)

            net_x, net_y, net_z = np.sum(bulk_displ_rel, axis = 0)
            data.attributes["netMSDx_bulk_com"] = np.mean(net_x**2) #/ len(bulk_displ_rel) #netMSD = ionic diffusivity
            data.attributes["netMSDy_bulk_com"] = np.mean(net_y**2) #/ len(bulk_displ_rel)
            data.attributes["netMSDz_bulk_com"] = np.mean(net_z**2) #/ len(bulk_displ_rel)
            data.attributes["netMSDav_bulk_com"] = (net_x**2)+(net_y**2)+(net_z**2) / len(bulk_displ_rel)
        
            # The follwoing code implement cros_displacement algo to get offdiagonal compoenets of MSD 
            #sum_of_pos_xy_n = ((np.sum(np.abs(bulk_displ_rel[:,0])) + np.abs(np.sum(bulk_displ_rel[:,1])))**2) / len(bulk_displ_rel)
            #sum_of_pos_xz_n = ((np.sum(np.abs(bulk_displ_rel[:,0])) + np.abs(np.sum(bulk_displ_rel[:,2])))**2) / len(bulk_displ_rel)
            #sum_of_pos_yz_n = ((np.sum(np.abs(bulk_displ_rel[:,1])) + np.abs(np.sum(bulk_displ_rel[:,2])))**2) / len(bulk_displ_rel)
        
            #data.attributes["netMSDxy_bulk_com"] = 0.5*(sum_of_pos_xy_n - net_x**2/len(np.abs(bulk_displ_rel)) - net_y**2/len(np.abs(bulk_displ_rel))) #np.mean(cd_xy**2) #np.mean(np.array(cd_xy)**2, axis=1) 
            #data.attributes["netMSDxz_bulk_com"] = 0.5*(sum_of_pos_xz_n - net_x**2/len(np.abs(bulk_displ_rel)) - net_z**2/len(np.abs(bulk_displ_rel))) #np.mean(cd_xz**2) #np.mean(np.array(cd_xz)**2, axis=1) 
            #data.attributes["netMSDyz_bulk_com"] = 0.5*(sum_of_pos_yz_n - net_y**2/len(np.abs(bulk_displ_rel)) - net_z**2/len(np.abs(bulk_displ_rel)))
        
            data.attributes["Timestep"] = time
     
        pipeline.modifiers.append(modify3)

        export_file(pipeline, file_2, "txt/attr",
        columns=["Timestep", "netMSDx_bulk_com", "netMSDy_bulk_com", "netMSDz_bulk_com","netMSDav_bulk_com"],
        multiple_frames=True)   
        

class smooth_ovito:

    def twist_gb(input_files = ["/home/schulze/simulations/manual-structures/jochen_18/00-minimize/01-heat-523/02-diffusion"],filename='diffusivity-attributes-rolling_mean.txt'):

        #input_files = ["/home/schulze/simulations/manual-structures/jochen_18/00-minimize/01-heat-523/02-diffusion"]

        # Chosing half-width of of grain boundary w in Angstrom
        # Check with Ovito!
        w = 10
        #start = timeit.timeit()        
        ###Import file

        for infile in input_files:

            list = glob.glob(f"{infile}/structure.dump.*")
            list = [file for file in list if int(file.split(".dump.")[-1])%1000 == 0]  #Only every 10th dump-file loaded --> 500 of 5000 dump files 	
            pipeline = import_file(list)
            n_frames = pipeline.source.num_frames
            Na_type_id = 5
            #print(infile)
            #timestep between frames in ps
            #ÜBERPRÜFEN OB ZEITABSTÄNDE STIMMEN
            dt = 0.002*1000    # 2 ps pro Frame --> 500 Frames eingeladen

            #Definitions for moving average
            half_frame = pipeline.source.num_frames//2
            every_nth_frame = 1

            # Set element

            def setup_particle_types(frame, data):
                types = data.particles_.particle_types_
                types.type_by_id_(1).name = "Zr"
                types.type_by_id_(2).name = "Si"    
                types.type_by_id_(3).name = "P"
                types.type_by_id_(4).name = "O"
                types.type_by_id_(4).radius = 0.3
                types.type_by_id_(5).name = "Na"
            pipeline.modifiers.append(setup_particle_types)


            # Set charges and mass

            def set_charge_mass(frame: int, data: DataCollection):

                charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
                mass_dict = {"Zr": 91.224,
                                "Si": 28.086,
                                "P": 30.974,
                                "O": 15.999,
                                "Na": 22.99}
                charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
                mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
                ptypes = data.particles_.particle_types_
                for i in range(data.particles.count):
                    charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                    mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]

            pipeline.modifiers.append(set_charge_mass)


            # Displacement vectors
            displ = CalculateDisplacementsModifier(minimum_image_convention = False, use_frame_offset = True )
            pipeline.modifiers.append(displ)


            # Calculate 'Relative Displacement' of host matrix

            def rel_displacement(frame, data, type_id_Na = 5):
                #with respect to center of mass of host matrix                                                                                                                                                                                        
                ptypes = data.particles['Particle Type']
                displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
                displ = data.particles['Displacement']
                masses = data.particles['Mass']
                #calculate center of mass displacement of host matrix                                                                                                                                                                                 
                total_mass_matrix = np.sum( masses[ (ptypes != type_id_Na) ] )
                sum = 0.0
                for i in range(len(displ_matrix)):
                        sum +=displ_matrix[i] * masses[i]
                com_displ_host_matrix = sum/total_mass_matrix
                data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
                data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
            import functools
            pipeline.modifiers.append(functools.partial(rel_displacement, type_id_Na = 5))

            ###Custom python function that calculates and returns the msd_raw, netmsd_raw, msd_com, netmsd_com of atoms with a specific particle type id

    
            ########## neuMS
            def define_gb_grain(frame: int, data: DataCollection, w = w, type_id_Na = 5):  ### "frame: int," benötigt?
                #Select GB: Choice has to be double-checked in Ovito!
                pos = data.particles.positions
                Selection = data.particles_.create_property("Selection")
                Selection[(pos[:,2] < data.cell[2][3] + w) ] = 1
                Selection[(pos[:,2] > data.cell[2][3] + data.cell[2][2]/2. - w) & (pos[:,2] < data.cell[2][3] + data.cell[2][2]/2. + w) ] = 1
                Selection[(pos[:,2] > data.cell[2][3] + data.cell[2][2] - w) ] = 1

                data.apply(FreezePropertyModifier(source_property = 'Selection', 
                                          destination_property = 'Select0',
                                          freeze_at = 0))

                ptype = data.particles.particle_type
                #only Na atoms
                Selection = data.particles["Select0"]
                #print(type(Selection))

            pipeline.modifiers.append(define_gb_grain)

    
    
            ########## neuMS
            def Calculate_MSD_gb_raw(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                Selection = data.particles["Select0"] #####neuMS
                gb_displ = data.particles['Displacement'][(ptypes == type_id_Na) & (Selection)]
                return np.mean(gb_displ**2, axis=0)

            def Calculate_netMSD_gb_raw(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                Selection = data.particles["Select0"] #####neuMS
                gb_displ = data.particles['Displacement'][(ptypes == type_id_Na) & (Selection)]
                return np.mean(gb_displ, axis = 0)**2

            def Calculate_MSD_gb_com(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                Selection = data.particles["Select0"] #####neuMS
                gb_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na) & (Selection)]
                return np.mean(gb_displ_rel**2, axis = 0)

            def Calculate_netMSD_gb_com(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                Selection = data.particles["Select0"] #####neuMS
                gb_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na) & (Selection)]
                return np.mean(gb_displ_rel, axis = 0)**2
    
            def Calculate_MSD_grain_raw(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                #Selection = data.particles["Select0"] #####neuMS
                grain_displ = data.particles['Displacement'][(ptypes == type_id_Na)]
                return np.mean(grain_displ**2, axis=0)

            def Calculate_netMSD_grain_raw(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                #Selection = data.particles["Select0"] #####neuMS
                grain_displ = data.particles['Displacement'][(ptypes == type_id_Na)]
                return np.mean(grain_displ, axis = 0)**2

            def Calculate_MSD_grain_com(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                #Selection = data.particles["Select0"] #####neuMS
                grain_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na)]
                return np.mean(grain_displ_rel**2, axis = 0)

            def Calculate_netMSD_grain_com(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                #Selection = data.particles["Select0"] #####neuMS
                grain_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na)]
                return np.mean(grain_displ_rel, axis = 0)**2

            # Smoothing
    
            Na_MSD_gb = [np.array((0.0, 0.0, 0.0))]
            net_Na_MSD_gb = [np.array((0.0, 0.0, 0.0))]
            Na_MSD_com_gb = [np.array((0.0, 0.0, 0.0))]
            net_Na_MSD_com_gb = [np.array((0.0, 0.0, 0.0))]
    
            Na_MSD_grain = [np.array((0.0, 0.0, 0.0))]
            net_Na_MSD_grain = [np.array((0.0, 0.0, 0.0))]
            Na_MSD_com_grain = [np.array((0.0, 0.0, 0.0))]
            net_Na_MSD_com_grain = [np.array((0.0, 0.0, 0.0))]


            for i in range(every_nth_frame, half_frame, every_nth_frame):
                print(i)
            
                #end = timeit.timeit()
                #print(f"Outer for loop_begin: {i}, elapsed time: {(end - start)/60.:.2f} min")
            
                displ.frame_offset = -i
            
                Na_msd_gb = np.array((0.0, 0.0, 0.0))
                net_Na_msd_gb = np.array((0.0, 0.0, 0.0))
                Na_msd_com_gb = np.array((0.0, 0.0, 0.0))
                net_Na_msd_com_gb = np.array((0.0, 0.0, 0.0))
            
                Na_msd_grain = np.array((0.0, 0.0, 0.0))
                net_Na_msd_grain = np.array((0.0, 0.0, 0.0))
                Na_msd_com_grain = np.array((0.0, 0.0, 0.0))
                net_Na_msd_com_grain = np.array((0.0, 0.0, 0.0))

                for frame in range(i, n_frames): 
                       data = pipeline.compute(frame)
                   
                       Na_msd_gb += Calculate_MSD_gb_raw(data, Na_type_id)
                       net_Na_msd_gb += Calculate_netMSD_gb_raw(data, Na_type_id)
                       Na_msd_com_gb += Calculate_MSD_gb_com(data, Na_type_id)
                       net_Na_msd_com_gb += Calculate_netMSD_gb_com(data, Na_type_id)
                   
                       Na_msd_grain += Calculate_MSD_grain_raw(data, Na_type_id)
                       net_Na_msd_grain += Calculate_netMSD_grain_raw(data, Na_type_id)
                       Na_msd_com_grain += Calculate_MSD_grain_com(data, Na_type_id)
                       net_Na_msd_com_grain += Calculate_netMSD_grain_com(data, Na_type_id)

            
                Na_MSD_gb.append(Na_msd_gb/(n_frames - i))
                net_Na_MSD_gb.append(net_Na_msd_gb /(n_frames -i ))
                Na_MSD_com_gb.append(Na_msd_com_gb/(n_frames - i))
                net_Na_MSD_com_gb.append(net_Na_msd_com_gb /(n_frames -i ))
            
                Na_MSD_grain.append(Na_msd_grain/(n_frames - i))
                net_Na_MSD_grain.append(net_Na_msd_grain /(n_frames -i ))
                Na_MSD_com_grain.append(Na_msd_com_grain/(n_frames - i))
                net_Na_MSD_com_grain.append(net_Na_msd_com_grain /(n_frames -i ))


            #Export results to txt file
            t = np.arange(len(Na_MSD_grain), dtype = float)
            t*=(dt*every_nth_frame)
            np.savetxt( f"{infile}/diffusivity-attributes-rolling_mean.txt", np.column_stack( (t, Na_MSD_gb[:][0], Na_MSD_gb[:][1], Na_MSD_gb[:][2], 
                                                                                     net_Na_MSD_gb[:][0], net_Na_MSD_gb[:][1], net_Na_MSD_gb[:][2], 
                                                                                     Na_MSD_com_gb[:][0], Na_MSD_com_gb[:][1], Na_MSD_com_gb[:][2], 
                                                                                     net_Na_MSD_com_gb[:][0], net_Na_MSD_com_gb[:][1], net_Na_MSD_com_gb[:][2],
                                                                                     Na_MSD_grain[:][0], Na_MSD_grain[:][1], Na_MSD_grain[:][2], 
                                                                                     net_Na_MSD_grain[:][0], net_Na_MSD_grain[:][1], net_Na_MSD_grain[:][2], 
                                                                                     Na_MSD_com_grain[:][0], Na_MSD_com_grain[:][1], Na_MSD_com_grain[:][2], 
                                                                                     net_Na_MSD_com_grain[:][0], net_Na_MSD_com_grain[:][1], net_Na_MSD_com_grain[:][2],
                                                                                    )), delimiter = " ", header = "Timestep(ps),MSDx_gb_raw(A^2),MSDy_gb_raw(A^2),MSDz_gb_raw(A^2),netMSDx_gb_raw(A^2),netMSDy_gb_raw(A^2),netMSDz_gb_raw(A^2),MSDx_gb_com(A^2),MSDy_gb_com(A^2),MSDz_gb_com(A^2),netMSDx_gb_com(A^2),netMSDy_gb_com(A^2),netMSDz_gb_com(A^2),MSDx_grain_raw(A^2),MSDy_grain_raw(A^2),MSDz_grain_raw(A^2),netMSDx_grain_raw(A^2),netMSDy_grain_raw(A^2),netMSDz_grain_raw(A^2),MSDx_grain_com(A^2),MSDy_grain_com(A^2),MSDz_grain_com(A^2),netMSDx_grain_com(A^2),netMSDy_grain_com(A^2),netMSDz_grain_com(A^2)")


    def tilt_gb(input_files = ["/nfshome/deshmukh/vaibhav/schulze_project/gb/IDp03_aimsgbCK/00-minimize/01-heat-523/02-diffusion"],nth_frame=100,filename='diffusivity-attributes-rolling_mean.txt'):

        #input_files = ["/home/schulze/simulations/manual-structures/jochen_18/00-minimize/01-heat-523/02-diffusion"]

        # Chosing half-width of of grain boundary w in Angstrom
        # Check with Ovito!
        w = 10
        #start = timeit.timeit()        
        ###Import file
        
        for infile in input_files:
            print(infile)
            print(infile + '/' + filename)
            #list = glob.glob(f"{infile}/structure.dump.*")
            #list = [file for file in list if int(file.split(".dump.")[-1])%1000 == 0]  #Only every 10th dump-file loaded --> 500 of 5000 dump files 	
            pipeline = import_file(f"{infile}/structure.dump.*")
            n_frames = pipeline.source.num_frames
            Na_type_id = 5
            #print(infile)
            #timestep between frames in ps
            #ÜBERPRÜFEN OB ZEITABSTÄNDE STIMMEN
            dt = 0.002*1000    # 2 ps pro Frame --> 500 Frames eingeladen

            #Definitions for moving average
            half_frame = pipeline.source.num_frames//2
            every_nth_frame = nth_frame

            # Set element

            def setup_particle_types(frame, data):
                types = data.particles_.particle_types_
                types.type_by_id_(1).name = "Zr"
                types.type_by_id_(2).name = "Si"    
                types.type_by_id_(3).name = "P"
                types.type_by_id_(4).name = "O"
                types.type_by_id_(4).radius = 0.3
                types.type_by_id_(5).name = "Na"
            pipeline.modifiers.append(setup_particle_types)


            # Set charges and mass

            def set_charge_mass(frame: int, data: DataCollection):

                charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
                mass_dict = {"Zr": 91.224,
                                "Si": 28.086,
                                "P": 30.974,
                                "O": 15.999,
                                "Na": 22.99}
                charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
                mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
                ptypes = data.particles_.particle_types_
                for i in range(data.particles.count):
                    charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
                    mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]

            pipeline.modifiers.append(set_charge_mass)


            # Displacement vectors
            displ = CalculateDisplacementsModifier(minimum_image_convention = False, use_frame_offset = True )
            pipeline.modifiers.append(displ)


            # Calculate 'Relative Displacement' of host matrix

            def rel_displacement(frame, data, type_id_Na = 5):
                #with respect to center of mass of host matrix                                                                                                                                                                                        
                ptypes = data.particles['Particle Type']
                displ_matrix = data.particles['Displacement'][ (ptypes != type_id_Na) ]
                displ = data.particles['Displacement']
                masses = data.particles['Mass']
                #calculate center of mass displacement of host matrix                                                                                                                                                                                 
                total_mass_matrix = np.sum( masses[ (ptypes != type_id_Na) ] )
                sum = 0.0
                for i in range(len(displ_matrix)):
                        sum +=displ_matrix[i] * masses[i]
                com_displ_host_matrix = sum/total_mass_matrix
                data.attributes['com_displ_host_matrix'] = com_displ_host_matrix
                data.particles_.create_property('Relative Displacement', data = displ - com_displ_host_matrix)
            import functools
            pipeline.modifiers.append(functools.partial(rel_displacement, type_id_Na = 5))

            ###Custom python function that calculates and returns the msd_raw, netmsd_raw, msd_com, netmsd_com of atoms with a specific particle type id

    
            ########## neuMS
            def define_gb_grain(frame: int, data: DataCollection, w = w, type_id_Na = 5):  ### "frame: int," benötigt?
                #Select GB: Choice has to be double-checked in Ovito!
                pos = data.particles.positions
                Selection = data.particles_.create_property("Selection")
                pos = data.particles.positions
                Selection[(pos[:,0] < data.cell[0][3] + w) ] = 1
                Selection[(pos[:,0] > data.cell[0][3] + data.cell[0][0]/2. - w) & (pos[:,0] < data.cell[0][3] + data.cell[0][0]/2. + w) ] = 1
                Selection[(pos[:,0] > data.cell[0][3] + data.cell[0][0] - w) ] = 1

                data.apply(FreezePropertyModifier(source_property = 'Selection', 
                                          destination_property = 'Select0',
                                          freeze_at = 0))

                ptype = data.particles.particle_type
                #only Na atoms
                Selection = data.particles["Select0"]
                #print(type(Selection))

            pipeline.modifiers.append(define_gb_grain)

    
    
            ########## neuMS
            def Calculate_MSD_gb_raw(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                Selection = data.particles["Select0"] #####neuMS
                gb_displ = data.particles['Displacement'][(ptypes == type_id_Na) & (Selection)]
                return np.mean(gb_displ**2, axis=0)

            def Calculate_netMSD_gb_raw(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                Selection = data.particles["Select0"] #####neuMS
                gb_displ = data.particles['Displacement'][(ptypes == type_id_Na) & (Selection)]
                return np.mean(gb_displ, axis = 0)**2

            def Calculate_MSD_gb_com(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                Selection = data.particles["Select0"] #####neuMS
                gb_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na) & (Selection)]
                return np.mean(gb_displ_rel**2, axis = 0)

            def Calculate_netMSD_gb_com(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                Selection = data.particles["Select0"] #####neuMS
                gb_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na) & (Selection)]
                return np.mean(gb_displ_rel, axis = 0)**2
    
            def Calculate_MSD_grain_raw(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                #Selection = data.particles["Select0"] #####neuMS
                grain_displ = data.particles['Displacement'][(ptypes == type_id_Na)]
                return np.mean(grain_displ**2, axis=0)

            def Calculate_netMSD_grain_raw(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                #Selection = data.particles["Select0"] #####neuMS
                grain_displ = data.particles['Displacement'][(ptypes == type_id_Na)]
                return np.mean(grain_displ, axis = 0)**2

            def Calculate_MSD_grain_com(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                #Selection = data.particles["Select0"] #####neuMS
                grain_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na)]
                return np.mean(grain_displ_rel**2, axis = 0)

            def Calculate_netMSD_grain_com(data, type_id_Na):
                ptypes = data.particles["Particle Type"]
                #Selection = data.particles["Select0"] #####neuMS
                grain_displ_rel = data.particles['Relative Displacement'][(ptypes == type_id_Na)]
                return np.mean(grain_displ_rel, axis = 0)**2

            # Smoothing
    
            Na_MSD_gb = [np.array((0.0, 0.0, 0.0))]
            net_Na_MSD_gb = [np.array((0.0, 0.0, 0.0))]
            Na_MSD_com_gb = [np.array((0.0, 0.0, 0.0))]
            net_Na_MSD_com_gb = [np.array((0.0, 0.0, 0.0))]
    
            Na_MSD_grain = [np.array((0.0, 0.0, 0.0))]
            net_Na_MSD_grain = [np.array((0.0, 0.0, 0.0))]
            Na_MSD_com_grain = [np.array((0.0, 0.0, 0.0))]
            net_Na_MSD_com_grain = [np.array((0.0, 0.0, 0.0))]


            for i in range(every_nth_frame, half_frame, every_nth_frame):
                print(i)
            
                #end = timeit.timeit()
                #print(f"Outer for loop_begin: {i}, elapsed time: {(end - start)/60.:.2f} min")
            
                displ.frame_offset = -i
            
                Na_msd_gb = np.array((0.0, 0.0, 0.0))
                net_Na_msd_gb = np.array((0.0, 0.0, 0.0))
                Na_msd_com_gb = np.array((0.0, 0.0, 0.0))
                net_Na_msd_com_gb = np.array((0.0, 0.0, 0.0))
            
                Na_msd_grain = np.array((0.0, 0.0, 0.0))
                net_Na_msd_grain = np.array((0.0, 0.0, 0.0))
                Na_msd_com_grain = np.array((0.0, 0.0, 0.0))
                net_Na_msd_com_grain = np.array((0.0, 0.0, 0.0))

                for frame in range(i, n_frames): 
                       data = pipeline.compute(frame)
                   
                       Na_msd_gb += Calculate_MSD_gb_raw(data, Na_type_id)
                       net_Na_msd_gb += Calculate_netMSD_gb_raw(data, Na_type_id)
                       Na_msd_com_gb += Calculate_MSD_gb_com(data, Na_type_id)
                       net_Na_msd_com_gb += Calculate_netMSD_gb_com(data, Na_type_id)
                   
                       Na_msd_grain += Calculate_MSD_grain_raw(data, Na_type_id)
                       net_Na_msd_grain += Calculate_netMSD_grain_raw(data, Na_type_id)
                       Na_msd_com_grain += Calculate_MSD_grain_com(data, Na_type_id)
                       net_Na_msd_com_grain += Calculate_netMSD_grain_com(data, Na_type_id)

            
                Na_MSD_gb.append(Na_msd_gb/(n_frames - i))
                net_Na_MSD_gb.append(net_Na_msd_gb /(n_frames -i ))
                Na_MSD_com_gb.append(Na_msd_com_gb/(n_frames - i))
                net_Na_MSD_com_gb.append(net_Na_msd_com_gb /(n_frames -i ))
            
                Na_MSD_grain.append(Na_msd_grain/(n_frames - i))
                net_Na_MSD_grain.append(net_Na_msd_grain /(n_frames -i ))
                Na_MSD_com_grain.append(Na_msd_com_grain/(n_frames - i))
                net_Na_MSD_com_grain.append(net_Na_msd_com_grain /(n_frames -i ))


            #Export results to txt file
            t = np.arange(len(Na_MSD_grain), dtype = float)
            t*=(dt*every_nth_frame)

            #np.savetxt(infile + '/' + filename, np.column_stack( (t, Na_msd_gb[0], Na_msd_gb[1], Na_msd_gb[2], net_Na_msd_gb[0], net_Na_msd_gb[1], net_Na_msd_gb[2], Na_msd_com_gb[0], Na_msd_com_gb[1], Na_msd_com_gb[2], net_Na_msd_com_gb[0], net_Na_msd_com_gb[1], net_Na_msd_com_gb[2], Na_msd_grain[0], Na_msd_grain[1], Na_msd_grain[2], net_Na_msd_grain[0], net_Na_msd_grain[1], net_Na_msd_grain[2], Na_msd_com_grain[0], Na_msd_com_grain[1], Na_msd_com_grain[2], net_Na_msd_com_grain[0], net_Na_msd_com_grain[1], net_Na_msd_com_grain[2] )), delimiter = " ")#, header = "Timestep(ps),MSDx_gb_raw(A^2),MSDy_gb_raw(A^2),MSDz_gb_raw(A^2),netMSDx_gb_raw(A^2),netMSDy_gb_raw(A^2),netMSDz_gb_raw(A^2),MSDx_gb_com(A^2),MSDy_gb_com(A^2),MSDz_gb_com(A^2),netMSDx_gb_com(A^2),netMSDy_gb_com(A^2),netMSDz_gb_com(A^2),MSDx_grain_raw(A^2),MSDy_grain_raw(A^2),MSDz_grain_raw(A^2),netMSDx_grain_raw(A^2),netMSDy_grain_raw(A^2),netMSDz_grain_raw(A^2),MSDx_grain_com(A^2),MSDy_grain_com(A^2),MSDz_grain_com(A^2),netMSDx_grain_com(A^2),netMSDy_grain_com(A^2),netMSDz_grain_com(A^2)")
            np.savetxt( infile + '/' + filename, np.column_stack( (t, Na_MSD_gb[:][0], Na_MSD_gb[:][1], Na_MSD_gb[:][2], 
                                                                                     net_Na_MSD_gb[:][0], net_Na_MSD_gb[:][1], net_Na_MSD_gb[:][2], 
                                                                                     Na_MSD_com_gb[:][0], Na_MSD_com_gb[:][1], Na_MSD_com_gb[:][2], 
                                                                                     net_Na_MSD_com_gb[:][0], net_Na_MSD_com_gb[:][1], net_Na_MSD_com_gb[:][2],
                                                                                     Na_MSD_grain[:][0], Na_MSD_grain[:][1], Na_MSD_grain[:][2], 
                                                                                     net_Na_MSD_grain[:][0], net_Na_MSD_grain[:][1], net_Na_MSD_grain[:][2], 
                                                                                     Na_MSD_com_grain[:][0], Na_MSD_com_grain[:][1], Na_MSD_com_grain[:][2], 
                                                                                     net_Na_MSD_com_grain[:][0], net_Na_MSD_com_grain[:][1], net_Na_MSD_com_grain[:][2],
                                                                                    )), delimiter = " ", header = "Timestep(ps),MSDx_gb_raw(A^2),MSDy_gb_raw(A^2),MSDz_gb_raw(A^2),netMSDx_gb_raw(A^2),netMSDy_gb_raw(A^2),netMSDz_gb_raw(A^2),MSDx_gb_com(A^2),MSDy_gb_com(A^2),MSDz_gb_com(A^2),netMSDx_gb_com(A^2),netMSDy_gb_com(A^2),netMSDz_gb_com(A^2),MSDx_grain_raw(A^2),MSDy_grain_raw(A^2),MSDz_grain_raw(A^2),netMSDx_grain_raw(A^2),netMSDy_grain_raw(A^2),netMSDz_grain_raw(A^2),MSDx_grain_com(A^2),MSDy_grain_com(A^2),MSDz_grain_com(A^2),netMSDx_grain_com(A^2),netMSDy_grain_com(A^2),netMSDz_grain_com(A^2)")
