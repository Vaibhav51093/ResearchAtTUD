from pyiron import Project
from pyiron import ase_to_pyiron, pyiron_to_ase
import sys  
sys.path.insert(0, '/nfshome/deshmukh/vaibhav/scripts')
import analysis_msd as ms

li = [573,623,673,723,773]
for i in li:
    ms.new_ovito.charge_diff(path_line='/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_2_2/glass/glass_hena_2/glass_6_%s/dump_nvt_prod.out'%i,file_2='msd_na_big_glass_ch_%sk_6.txt'%i,step=2.00, dump_write=200, step_ovito=1) #right one
    
    