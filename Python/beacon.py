

    from gpatom.beacon import BEACON
    from gpatom.beacon.str_gen import RandomBranch
    from ase.calculators.emt import EMT
    from ase import Atoms
    from ase import Atoms
    from ase.io import read, write
    from gpaw import GPAW, PW
    import numpy as np
    import random
    
    # Assuming calc, calcparams, and other necessary imports are defined

    class SiPSwapGenerator:
        def __init__(self, base_structure):
            self.base_structure = base_structure

        def get(self):
            # Make a copy of the structure to keep the original unchanged
            new_structure = self.base_structure.copy()

            # Get indices of Si and P atoms
            si_indices = [atom.index for atom in new_structure if atom.symbol == 'Si']
            p_indices = [atom.index for atom in new_structure if atom.symbol == 'P']

            # Combine and shuffle indices
            combined_indices = si_indices + p_indices
            random.shuffle(combined_indices)

            # Reassign atom types
            for i, idx in enumerate(combined_indices):
                if i < len(si_indices):
                    new_structure[idx].symbol = 'Si'
                else:
                    new_structure[idx].symbol = 'P'

            return new_structure


    # Load your initial structure (replace 'initial_structure_file' with your file path)
    initial_structure = read('/nfshome/deshmukh/vaibhav/NaSICON_dft/database/pbe_structure_10_hdf5/pbe_structure_10/CONTCAR', format='vasp')

    # Initialize your custom structure generator
    sgen = SiPSwapGenerator(initial_structure)

    calc = GPAW
    calc_params = dict(mode=PW(300),
                       xc='PBE',
                       symmetry='off',
                       kpts={'size': (2, 2, 2)},
                       txt='gpaw.txt')

    #calc = EMT
    #calc_params = dict()

    # Define the GEOFEE optimizer with your custom structure generator
    go = BEACON(calc=calc,
                calcparams=calc_params,
                sgen=sgen,
                nsur=1,            # Number of surface atoms (adjust as necessary)
                nbest=5,           # Number of best structures to keep
                nrattle=0,         # Number of structures to rattle
                ndft=40, 
                rng = np.random.RandomState(90))           # Number of DFT calculations to perform

    # Run the optimization
    go.run()
    
