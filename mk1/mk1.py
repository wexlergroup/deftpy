import pymatgen.core.structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view
from pymatgen.io.cif import CifParser

# from mp_api.client import MPRester
# with MPRester(api_key="C4VoKwbKXbVC9dk85Uq5G5tig1LWbPKu") as mpr:
#     data = mpr.materials.get_data_by_id("mp-4019")

class crystalClass:
    def file_readin(filething):
        if filething.startswith('mp-'):
            with MPRester(api_key="something") as mpr:
                data = mpr.materials.get_data_by_id(filething)

        else:
            structure = pymatgen.core.Structure.from_file(filething)
            visualized = AseAtomsAdaptor.get_atoms(structure)
            # Supported formats include CIF, POSCAR/CONTCAR, CHGCAR, LOCPOT, 
            # vasprun.xml, CSSR, Netcdf and pymatgen’s JSON-serialized structures.
        return(view(visualized))
    file_readin('mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt')



