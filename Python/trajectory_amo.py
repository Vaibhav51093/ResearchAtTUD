from pyiron import Project
from pyiron import ase_to_pyiron, pyiron_to_ase
import sys  
sys.path.insert(0, '/nfshome/deshmukh/vaibhav/scripts')

from analysis_msd_2 import oviot


oviot.trajectory(path_line="/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/different_quench_2/boilot_443K_2_05_amo_s_11_rate_1_k_ps_hdf5/boilot_443K_2_05_amo_s_11_rate_1_k_ps/dump.out", outfile='last_struct_amo.xyz', frame=100)
