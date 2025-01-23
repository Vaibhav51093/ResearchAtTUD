# run in parallel
import ase.io
from gpaw import GPAW, PW
#
atoms = ase.io.read("POSCAR-001")
#change the volume
calc  = GPAW(mode=PW(200), kpts=(2, 2, 2), xc="PBE", txt="displacement-1.txt")
atoms.set_calculator(calc)
atoms.get_forces()
