from pyiron import Project, ase_to_pyiron
import matplotlib.pyplot as plt
import numpy as np
from pyiron import Project
from ase.io import read, write
from pyiron import ase_to_pyiron
import ase
import os
import time


# 1. Job transfer form remote 
class basic_operations:
    def job_transfer(path="/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_proto_glass",job_name='boilot_623K_2_05_cryt_random_523_k_ps_new'):
        pr = Project(path) 
        job= pr[job_name]
        job.transfer_from_remote()
        job.compress()

    # 2. Decompress file     
    def job_d_comp(path="/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_proto_glass",job_name='boilot_623K_2_05_cryt_random_523_k_ps_new'):
        pr = Project(path) 
        job= pr[job_name]
        job.decompress()

    # 3. Compress file 
    def job_comp(path="/nfshome/deshmukh/pyiron/projects/NASICON/project/padone_pot/new_proto_glass",job_name='boilot_623K_2_05_cryt_random_523_k_ps_new'):
        pr = Project(path) 
        job= pr[job_name]
        job.compress()