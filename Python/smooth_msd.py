from pyiron import Project
from pyiron import ase_to_pyiron, pyiron_to_ase
import sys  
sys.path.insert(0, '/nfshome/deshmukh/vaibhav/scripts')
import analysis_smooth_msd as ms 
import analysis_msd as msd 

msd.smooth_ovito.tilt_gb(input_files = ["/nfshome/deshmukh/vaibhav/schulze_project/gb/IDp03_aimsgbCK/00-minimize/01-heat-523/02-diffusion"],nth_frame=1,filename='diffusivity-attributes-rolling_mean.txt')