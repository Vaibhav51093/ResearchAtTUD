from pyiron import Project
from pyiron import ase_to_pyiron, pyiron_to_ase
import sys  
sys.path.insert(0, '/nfshome/deshmukh/vaibhav/scripts')
import analysis_msd as ms

list = [43] #, 47, 4, 13, 42]
list_2 = [523, 573, 623, 673, 723, 773]

#for i in list:
for j in list_2:
#        #job = pr['hena_1_struct_eq_%sk_%s'%(i,j)]
    ms.new_ovito.charge_diff(path_line='/nfshome/deshmukh/pyiron/projects/NASICON/project/hena/hena_1/glass/hena_1_glass/glass_43_%sk/dump_nvt_prod.out'%(j),file_2='msd_na_big_glass_ch_%sk_43.txt'%(j),step=2.00, dump_write=200, step_ovito=1) #right one
    
    