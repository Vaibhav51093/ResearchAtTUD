from pyiron import Project
from pyiron import ase_to_pyiron, pyiron_to_ase
import sys  
sys.path.insert(0, '/nfshome/deshmukh/vaibhav/scripts')

from analysis_msd_2 import oviot


oviot.trajectory(path_line="/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/glass/hena_1_glass/glass_43_623k/dump_nvt_prod.out", outfile='last_struct_hena_1_amo.xyz', frame=100)
