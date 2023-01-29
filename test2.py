import pymatgen.core.structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view
from pymatgen.ext.matproj import MPRester

def blah(*args):
    structure = pymatgen.core.IStructure.from_file(args[0])
    get_neighbors = structure.get_all_neighbors_py()


blah("mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt")