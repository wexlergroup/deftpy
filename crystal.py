# standard library imports
from typing import Union, Any

# third-party imports
import numpy as np
import pandas as pd
from ase.visualize import view
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.structure_matcher import StructureMatcher as structure_matcher
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor

# local imports
import pymatgen.core.structure


class Crystal:
    """
    A class representing a crystal structure.

    Attributes:
        structure (Structure): The pymatgen Structure object representing the crystal structure.
        unique_oxygens (dict): A dictionary of unique oxygen sites in the crystal structure. 
        nonO (dict): A dictionary of non-oxygen sites in the crystal structure. 
                     
    Methods:
        __init__(self, structure: Structure): Initializes a new Crystal object with the given pymatgen Structure object.
        visualize(self): Opens a visualization window of the crystal structure using ASE's view function.
        unique_oxygen(self, structure: Structure): Identifies and stores the unique oxygen sites in the crystal 
                                                     structure as a dictionary of pymatgen Structure objects.
        non_oxygen(self, structure: Structure): Identifies and stores the non-oxygen sites in the crystal structure 
                                                 as a dictionary of pymatgen Site objects.
    """

    def __init__(self, structure: Structure):
        self.structure = structure

    def visualize(self):
        atoms = AseAtomsAdaptor.get_atoms(self.structure)
        view(atoms)

    def unique_oxygen(self, structure):
        structure.add_oxidation_state_by_guess()
        iO = []
        vac_structures = {}
        unique_oxygens = {}

        for i, site in enumerate(structure.sites):
            if site.specie.symbol == 'O':
                iO.append(i)

        for i in iO:
            new_struct = structure.copy()
            new_struct.remove_sites([i])
            vac_structures[i] = new_struct

        for i, vac_structure in vac_structures.items():
            if i == 0:
                unique_oxygens[i] = vac_structure
            n_dupl = 0
            for j, unique_structure in unique_oxygens.items():
                if structure_matcher.fit(vac_structure, unique_structure) == True:
                    n_dupl += 1
            if n_dupl == 0:
                unique_oxygens[i] = vac_structure
        
        self.unique_oxygens = unique_oxygens

    def non_oxygen(self, structure):
        nonO = {}
        for i, site in enumerate(structure.sites):
            if site.specie.symbol != 'O':
                nonO[i] = site
        self.nonO = nonO




def file_readin(filepath: str, visualized: bool = False) -> Union[Crystal, Any]:
    """
    Reads a crystal structure from a file and returns a `Crystal` object.

    Args:
        filepath (str): The path to the file containing the crystal structure.
        visualized (bool, optional): Whether to visualize the crystal structure using ASE. Defaults to False.

    Returns:
        Crystal or Any: The `Crystal` object representing the crystal structure, or `None` if the file cannot be read.
    """
    structure = pymatgen.core.Structure.from_file(filepath)
    crystal = Crystal(structure)
    if visualized:
        crystal.visualize()

    return(crystal)


def oxygen_cn(crystal: Crystal) -> dict:
    """
    Computes the coordination numbers and charges of the unique oxygen atoms in a crystal structure.
    
    Args:
        crystal (Crystal): The Crystal object to analyze.
    
    Returns:
        dict: A dictionary where the keys are the indices of the unique oxygen atoms in the
        structure, and the values are dictionaries containing information on the coordination
        numbers and charges of the atom's nearest neighbors.
    """
      
    structure = crystal.structure
    unique_oxygens = crystal.unique_oxygens
    
    struc_NN = {}
    for i in range(len(unique_oxygens.keys())):
         pos_arg = list(unique_oxygens.keys())[i]

    struc_NN[i] = CrystalNN(weighted_cn=True, cation_anion=True, porous_adjustment=False,
                                distance_cutoffs=None, x_diff_weight=0).get_nn_info(structure, pos_arg)

    return(struc_NN)

def non_oxygen_oxi_state(crystal: Crystal) -> dict:
    """
    Computes the oxidation states of the non-oxygen atoms in a crystal structure.
    
    Args:
        crystal (Crystal): The Crystal object to analyze.
        
    Returns:
        dict: A dictionary where the keys are the indices of the non-oxygen atoms in the
        structure, and the values are the oxidation states of the atoms.
    """
    structure = crystal.structure
    nonO = crystal.nonO
    nonO_oxi_state = {}
    for i in nonO:
         nonO_oxi_state[i] = structure[i].specie.add_charges_from_oxi_state_guesses
    return(nonO_oxi_state)


def crystal_data(crystal: Crystal):
    structure = crystal.structure
    unique_oxygens = crystal.unique_oxygens
    nonO = crystal.nonO

    data = []



        
    return('hi')
