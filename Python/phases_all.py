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

from ase.spacegroup import get_spacegroup
from matplotlib import pyplot as plt
plt.style.use(['vaibhz-sci','no-latex','high-vis','vibrant'])

from ase.io import read
from ase.db import connect
from ase.io.vasp import read_vasp, write_vasp

from ase.io import read
from ase.db import connect
from pymatgen.io.ase import AseAtomsAdaptor
import json
from pymatgen.ext.matproj import MPRester, Composition, Element
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry, GrandPotentialPhaseDiagram, GrandPotPDEntry, CompoundPhaseDiagram
#from pynter import SETTINGS
#from pynter.tools.utils import save_object_as_json, get_object_from_json
import requests
from ase.io import read, write
import re
import numpy as np
from ase.formula import Formula
import itertools
from sympy import symbols, Eq, solve, Rational
from ase.formula import Formula