import pymatgen.core.structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher
from ase.visualize import view
from pymatgen.ext.matproj import MPRester
#import CrystalNN
from pymatgen.analysis.local_env import CrystalNN
import pandas as pd
from pymatgen.core.composition import Composition


structure = pymatgen.core.Structure.from_file('mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt')
structure = pymatgen.core.Structure.from_file(args[0])
structure_matcher = pymatgen.analysis.structure_matcher.StructureMatcher()
structure.add_oxidation_state_by_guess()