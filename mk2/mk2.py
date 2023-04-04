import pymatgen.core.structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.core.composition import Composition

from ase.visualize import view

import pandas as pd
import numpy as np

class Crystal:

    def crystal_file(self, filepath):
        self.crystal_file = pymatgen.core.Structure.from_file(filepath)

    def unique_oxygens(self):
        '''returns a dictionary where keys are indices of
        unique oxygens in the structure and values are '''
        iO = []
        vac_structures = {}
        unique_structures = {}

        structure = self.crystal_file
        structure_matcher = pymatgen.analysis.structure_matcher.StructureMatcher()
        structure.add_oxidation_state_by_guess()

        for i, site in enumerate(structure.sites):
            if site.specie.symbol == 'O':
                iO.append(i)
        
        for i in iO:
            temp_struct = structure.copy()
            temp_struct.remove_sites([i])
            vac_structures[i] = temp_struct

        for i, vac_structure in vac_structures.items():
            if i == 0:
                unique_structures[i] = vac_structure
            n_dupl = 0
            for j, unique_structure in unique_structures.items():
                if structure_matcher.fit(vac_structure, unique_structure) == True:
                    n_dupl += 1
            if n_dupl == 0:
                unique_structures[i] = vac_structure
        return unique_structures

    def non_o_oxidation_states(self):
        '''returns a dictionary where'''
        notO = []
        oxi_states = {}
        structure = self.crystal_file
        
        for i, site in enumerate(structure.sites):
            if site.specie.symbol != 'O':
                notO.append(i)
        for j in notO:
            oxi_states[j] = structure[j].specie.add_charges_from_oxi_state_guesses

            return oxi_states
         
    
    # def unique_o_neighbors(self, filepath):
    #     return 
    
