import pymatgen.core.structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view
from pymatgen.io.cif import CifParser



blah = pymatgen.core.Structure.from_file('mk1/testPOSCAR.txt')
# Supported formats include CIF, POSCAR/CONTCAR, CHGCAR, LOCPOT, 
# vasprun.xml, CSSR, Netcdf and pymatgen’s JSON-serialized structures.
print(blah)



hi = AseAtomsAdaptor.get_atoms(blah)

view(hi)




