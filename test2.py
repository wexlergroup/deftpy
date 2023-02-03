import pymatgen.core.structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view
from pymatgen.ext.matproj import MPRester
import pymatgen.analysis.structure_matcher

structure = pymatgen.core.Structure.from_file("mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt")
iO = []
vac_structures = {}
for i, site in enumerate(structure.sites):
    if site.specie.symbol == 'O':
        iO.append(i)        

# print(iO)
# print('HIHIHIHHIHIHIHHIHIHIHHIHIHIHHIHIHIHHIHIHIHHIHIHIH')

for i in iO:
    new_struct = structure.copy()
    new_struct.remove_sites([i])
    vac_structures[i] = new_struct
# print(new_struct)

structure_matcher = pymatgen.analysis.structure_matcher.StructureMatcher()
for i in vac_structures:
    print(i)
    print(structure_matcher.fit(vac_structures[i], structure))
    print('HIHIHIHHIHIHIHHIHIHIHHIHIHIHHIHIHIHHIHIHIHHIHIHIH')
    


