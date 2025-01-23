# check the counterparts for each glass composition
import json
from pymatgen.ext.matproj import MPRester, Element
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, GrandPotentialPhaseDiagram, GrandPotPDEntry, CompoundPhaseDiagram
#from pynter import SETTINGS
#from pynter.tools.utils import save_object_as_json, get_object_from_json
import requests
from pymatgen.core.composition import Composition
from pymatgen.io.vasp.outputs import Outcar
from ase.io import read, write
from ase.io.vasp import read_vasp, write_vasp 
from pymatgen.ext.matproj import MPRester, Composition 
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PatchedPhaseDiagram, CompoundPhaseDiagram, ReactionDiagram
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.reaction_calculator import ComputedReaction

# Direct lammps script implemetation in pyiron 
from tempfile import tempdir
from pyiron import Project
import numpy as np
import pandas
from jinja2 import Template
import matplotlib.pyplot as plt 
import scipy.constants as sc
from scipy.integrate import cumtrapz
import os
from pyiron import ase_to_pyiron, pyiron_to_ase
from ase.io import read, write
from ase.formula import Formula

from ase.spacegroup import get_spacegroup
from matplotlib import pyplot as plt
from ase.db import connect
plt.style.use(['vaibhz-sci','no-latex','high-vis','vibrant'])

pr = Project("/nfshome/deshmukh/pyiron/projects/NASICON/project/frank/potetial_phase_diagram")

# Customizing the potentials (New potentials) [Adjust as per the Requrnments]
import pandas
padone_potential = pandas.DataFrame({
  'Name': ['NASICON_Pedone'],
  'Filename': [[]],
  'Model': ['Custom'],
  'Species': [['O', 'Na', 'Zr', 'Si', 'P']],
  'Config': [['atom_style full\n',      # Given function contains morse parametrs, i.e. it has bonds 
              '## create groups ###\n',
              'group O type 1\n',
              'group Na type 2\n',
              'group Zr type 3\n',
              'group Si type 4\n',
              'group P type 5\n',
              '\n',
              '## set charges - beside manually ###\n',
              'set group O charge -1.2000\n',
              'set group Na charge 0.6000\n',
              'set group Zr charge 2.4000\n',
              'set group Si charge 2.4000\n',
              'set group P charge 3.0000\n',
              '\n',
              'pair_style hybrid/overlay morse 15.0 mie/cut 15.0 coul/long 15.0 beck 15.0\n',
              'pair_coeff * * coul/long\n',
              'pair_coeff 1 2 beck 5.0 0 0 0 0\n',
              'pair_coeff 1 3 beck 1.0 0 0 0 0\n',
              'pair_coeff 1 4 beck 1.0 0 0 0 0\n',
              'pair_coeff 1 5 beck 1.0 0 0 0 0\n',
              'pair_coeff 1 1 beck 22.0 0 0 0 0\n',
              'pair_coeff 1 2 mie/cut 5.0 1.0 12.0 0\n',
              'pair_coeff 1 3 mie/cut 1.0 1.0 12.0 0\n',
              'pair_coeff 1 4 mie/cut 1.0 1.0 12.0 0\n',
              'pair_coeff 1 5 mie/cut 1.0 1.0 12.0 0\n',
              'pair_coeff 1 1 mie/cut 22.0 1.0 12.0 0\n',
              'pair_coeff 1 2 morse 0.023363 1.763867 3.006315\n',
              'pair_coeff 1 3 morse 0.206237 2.479675 2.436997\n',
              'pair_coeff 1 4 morse 0.340554 2.006700 2.100000\n',
              'pair_coeff 1 5 morse 0.831326 2.585833 1.800790\n',
              'pair_coeff 1 1 morse 0.042395 1.379316 3.618701\n',
              'kspace_style ewald 1.0e-8\n']]
})


db_2 = connect('/nfshome/deshmukh/vaibhav/share/new_data_jochen/NaSICON.ZrVacancies.equilibrated+relaxed.db')

for k in range(len(db_2)):
    print('-----------------------------------------------------------------------')
    #x = db_2.get(id=i+1).toatoms()
    #x.get_chemical_formula(empirical=False)
    if k+1 >= 257:
      bulk_struct = ase_to_pyiron(db_2.get(id=k+1).toatoms())
      atoms = pyiron_to_ase(bulk_struct)
      # check for charge neutrality if it 
      #print('Charge = ', 0.6*list(Formula(atoms.get_chemical_formula()).count().values())[0] - 1.2*list(Formula(atoms.get_chemical_formula()).count().values())[1] + 3*list(Formula(atoms.get_chemical_formula()).count().values())[2] + 2.4*list(Formula(atoms.get_chemical_formula()).count().values())[3] + 2.4*list(Formula(atoms.get_chemical_formula()).count().values())[4])
      job_name = 'phases_pot_struct_jochen_all_%s'%k
      job_lmp = pr.create.job.Lammps(job_name, delete_existing_job=True)
      job_lmp.structure = bulk_struct
      job_lmp.potential = padone_potential
      job_lmp.calc_minimize(pressure=0.0)
      job_lmp.server.queue = "normal_l2"    ### See the queues.yaml file for queue names. Standard queue on Lichtenberg is normal with a runtime limit of 24 hours.
      job_lmp.server.cores = 24             ### The number of cores you want for your job
      job_lmp.run()
    else:
      continue 
    

