import pymatgen.core.structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view
from pymatgen.ext.matproj import MPRester

# for using data from the Materials Project, enter a tuple with the 
# mp-id as the first element and your api key as the second element
def file_readin(*args):
    if args[0].startswith('mp-'):
        with MPRester(api_key=(args[1])) as mpr:
            data = mpr.materials.get_data_by_id(args[0]) 
    else:
        structure = pymatgen.core.Structure.from_file(args[0])
        visualized = AseAtomsAdaptor.get_atoms(structure)
        #works with cif, vasp, and poscar files

    return(view(visualized))
file_readin("mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt")

def blah(*args):
    structure = pymatgen.core.Structure.from_file(args[0])
    for x in structure.sites():
        print(x)
blah("mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt")


#test against vesta output


