from ase.io import read, write
from ase.optimize.bfgs import BFGS
from ase.calculators.vasp import Vasp


dft_input = read('/nfshome/deshmukh/vaibhav/frahunhofer_katja/nfep/pnma_NaFePO4_POSCAR.vasp', format='vasp')

calc = Vasp(prec='High',
            xc='PBE',
            lreal=False,
            kpts=k_p[0],
            encut=700,
            ismear=0,
            sigma=0.05,
            algo=48,
            nelm=200,
            nsw=200,
            ediff=1e-6,
            ediffg=-0.01,
            isif=3,
            atoms=dft_input)

opt = BFGS(calc, logfile='log.txt')
opt.run(fmax=0.01)
dft_input.get_potential_energy()
