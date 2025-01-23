import lammps_logfile 
#/home/vd80naku/miniconda3/lib/python3.9/site-packages/pandas/core/arrays/masked.py:60: UserWarning: Pandas requires version '1.3.6' or newer of 'bottleneck' (version '1.3.5' currently installed).
#from pandas.core import (
import numpy as np 
log_5_1 = lammps_logfile.File("log.lammps")
ly = log_5_1.get("Lz") 
print(np.mean(ly[-100:-1]))
