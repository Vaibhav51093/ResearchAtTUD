#---------------------------------------------------------------------------------#
# Official work of Vaibhav Deshmukh (TU-Darmstadt) 
# The script calculate the PDOS of the system and plot the PDOS
# It requires lammps dump file as input
# ASE package is required to run the script
#----------------------------------------------------------------------------#
from ase.io import read
import numpy as np 
import matplotlib.pyplot as plt
from scipy import fftpack
from scipy import signal
import numpy as np
from scipy.fftpack import fft, fftfreq
from ase.io.trajectory import Trajectory
from scipy.ndimage.filters import  gaussian_filter1d as gaussian
import time
from matplotlib import pyplot as plt
plt.style.use(['vaibhz-sci','no-latex'])
#---------------------------------------------------------------------------------#
# Record the start time
start_time = time.time()
# Print nice quotes about patience and perseverance
print('''---------------------------------------------------------------------------------------------------------------------''')
print('''Your calculations are in progress. Thanks for your patience and perseverance.''')

#---------------------------------------------------------------------------------#

# Main calculations and plotting
list_file = ['/nfshome/deshmukh/vaibhav/shanka_cluster/deshmukh2/nasi_2_random/523/dump_nvt_prod_3.out',
             '/nfshome/deshmukh/vaibhav/shanka_cluster/deshmukh2/nasi_2_random/573/dump_nvt_prod_3.out',
             '/nfshome/deshmukh/vaibhav/shanka_cluster/deshmukh2/nasi_2_random/623/dump_nvt_prod_3.out',
             '/nfshome/deshmukh/vaibhav/shanka_cluster/deshmukh2/nasi_2_random/673/dump_nvt_prod_3.out',
             '/nfshome/deshmukh/vaibhav/shanka_cluster/deshmukh2/nasi_2_random/723/dump_nvt_prod_3.out',
             '/nfshome/deshmukh/vaibhav/shanka_cluster/deshmukh2/nasi_2_random/773/dump_nvt_prod_3.out',]

# Timesteps
time_step = 2*200 # fs

# File names for saving the PDOS
total_file_name = 'pdos_overall_nasi_0.png'
na_file_name = 'pdos_na_nasi_0.png'

# Text file names for saving the max freq
text_file_t = 'max_freq_overall_nasi_0.txt'
text_file_na = 'max_freq_na_nasi_0.txt'

#---------------------------------------------------------------------------------#

def pdos(V, dt):
    """
    Calculate the phonon density of states from a trajectory of
    velocities (power spectrum of the velocity auto-correlation
    function).

    Parameters
    ----------
    V : :obj:`numpy.ndarray`
        (dims N x T) velocities of N degrees of freedom for
        trajetory of length T
    dt : float
        time between steps in trajectory (fs)

    Returns
    -------
    freq : :obj:`numpy.ndarray`
        (dims T) Frequencies (cm^-1)
    pdos : :obj:`numpy.ndarray`
        (dims T) Density of states (a.u.)
    """

    n_steps = V.shape[1]

    # mean velocity auto-correlation for all degrees of freedom
    vac2 = [np.correlate(v, v, 'full') for v in V]
    vac2 /= np.linalg.norm(vac2, axis=1)[:, None]
    vac2 = np.mean(vac2, axis=0)

    # power spectrum (phonon density of states)
    pdos = np.abs(fft(vac2))**2
    pdos /= np.linalg.norm(pdos) / 2 # spectrum is symmetric

    freq = fftfreq(2*n_steps-1, dt) * 33356.4095198152 * 0.0299792458 # cm^-1 to THz

    return freq[:n_steps], pdos[:n_steps]


def pdos_from_trajectory(file='/nfshome/deshmukh/vaibhav/shanka_cluster/deshmukh2/nas_0/523/dump_nvt_prod_3.out', dt=2*200):
    """
    Calculate the phonon density of states from a trajectory of
    velocities (power spectrum of the velocity auto-correlation
    function).

    Parameters
    ----------
    traj : :obj:`ase.io.trajectory.Trajectory`
        trajectory of atoms
    dt : float
        time between steps in trajectory (fs)

    Returns
    -------
    freq : :obj:`numpy.ndarray`
        (dims T) Frequencies (cm^-1)
    pdos : :obj:`numpy.ndarray`
        (dims T) Density of states (a.u.)
    """
    traj = read(file, index=':', format='lammps-dump-text')         # Read the trajectory file
    v_nasi = np.array([atoms.get_velocities() for atoms in traj])   # Get the velocities of all the atoms
    na_steps = v_nasi.shape[0]                                      # Number of steps in the trajectory
    v_nasi = v_nasi.reshape(na_steps, -1).T                         # Reshape the velocities
    
    freq, pdos_t = pdos(v_nasi, dt)                                 # Calculate the PDOS   

    v_na = []
    for atom in traj:
        na_indices = [i for i, atom in enumerate(atom) if atom.symbol == 'He']                  # He = Na :- The notations are not as per my stuff 
        vel_x = atom.get_velocities()[na_indices]
        v_na.append(vel_x)
    # Reshape the velocities    
    v_na = np.array(v_na)
    v_na = v_na.reshape(na_steps, -1).T
    freq_na, pdos_na = pdos(v_na, dt)
    
    return freq, gaussian(pdos_t, sigma=50), freq_na, gaussian(pdos_na, sigma=50)


# function to plot the PDOS for overall system and Na atoms

def plot_pdos(freq, pdos, file_name, label):
    
    plt.rcParams['figure.figsize'] = [6, 6]
    plt.plot(freq, pdos, label='%s'%label)
    plt.xlabel('Frequency (THz)')
    plt.ylabel('Density of states (a.u.)')
    plt.legend()
    plt.savefig(file_name)
    plt.show()
    plt.close() 
    
def plot_pdos_na(freq, pdos, file_name, label):
    
    plt.rcParams['figure.figsize'] = [6, 6]
    plt.plot(freq, pdos, label='%s'%label)
    plt.xlabel('Frequency (THz)')
    plt.ylabel('Density of states (a.u.)')
    plt.legend()
    plt.savefig(file_name)
    plt.show()
    plt.close()
    
    
#---------------------------------------------------------------------------------#

temp = [523, 573, 623, 673, 723, 773]
fr_1 = []
fr_na = []
pdos_1 = []
pdos_na = []
for file in list_file:
    freq, pdos_t, freq_na, pdos_n = pdos_from_trajectory(file, dt=time_step)
    fr_1.append(freq)
    fr_na.append(freq_na)
    pdos_1.append(pdos_t)
    pdos_na.append(pdos_n)
    
# Plotting the PDOS for overall system and Na atoms

# Plotting the PDOS for overall system and Na atoms

plt.rcParams['figure.figsize'] = [6, 6]
plt.plot(fr_1[0], pdos_1[0], label='T = %s K'%temp[0])
plt.plot(fr_1[1], pdos_1[1], label='T = %s K'%temp[1])
plt.plot(fr_1[2], pdos_1[2], label='T = %s K'%temp[2])
plt.plot(fr_1[3], pdos_1[3], label='T = %s K'%temp[3])
plt.plot(fr_1[4], pdos_1[4], label='T = %s K'%temp[4])
plt.plot(fr_1[5], pdos_1[5], label='T = %s K'%temp[5])
plt.xlabel('Frequency (THz)')
plt.ylabel('Density of states (a.u.)')
plt.grid(True, ls="--")
plt.legend()
plt.savefig(total_file_name)
plt.show()
plt.close()

# Plotting the PDOS for overall system and Na atoms

plt.rcParams['figure.figsize'] = [6, 6]
plt.plot(fr_na[0], pdos_na[0], label='T = %s K'%temp[0])
plt.plot(fr_na[1], pdos_na[1], label='T = %s K'%temp[1])
plt.plot(fr_na[2], pdos_na[2], label='T = %s K'%temp[2])
plt.plot(fr_na[3], pdos_na[3], label='T = %s K'%temp[3])
plt.plot(fr_na[4], pdos_na[4], label='T = %s K'%temp[4])
plt.plot(fr_na[5], pdos_na[5], label='T = %s K'%temp[5])
plt.xlabel('Frequency (THz)')
plt.ylabel('Density of states (a.u.)')
plt.grid(True, ls="--")
plt.legend()
plt.savefig(na_file_name)
plt.show()
plt.close()
# ---------------------------------------------------------------------------------#
# Write text file showing temp vs max freq in THz

with open(text_file_t, 'w') as f:
    for i in range(len(temp)):
        f.write('%s %s\n'%(temp[i], fr_1[i][np.argmax(pdos_1[i])]))

with open(text_file_na, 'w') as f:
    for i in range(len(temp)):
        f.write('%s %s\n'%(temp[i], fr_na[i][np.argmax(pdos_na[i])]))

#---------------------------------------------------------------------------------#
# Record the end time
end_time = time.time()
# Calculate the elapsed time
elapsed_time = end_time - start_time
# Priting the nice quotes about patience and perseverance
print('''Your calculations are completed. It took %.2f seconds.'''%elapsed_time)    
print('''---------------------------------------------------------------------------------------------------------------------''')
