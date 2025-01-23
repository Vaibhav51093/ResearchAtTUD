from pyiron import Project
from pyiron import ase_to_pyiron, pyiron_to_ase
import sys  
sys.path.insert(0, '/nfshome/deshmukh/vaibhav/scripts')

from analysis_msd_2 import oviot


oviot.trajectory(path_line="/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/minimization/hena_1_struct_eq_big_623k_43_hdf5/hena_1_struct_eq_big_623k_43/dump_nvt_prod.out", outfile='last_struct_hena_1_cryt.xyz', frame=100)
