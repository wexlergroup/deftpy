from pymatgen.io.cif import CifParser
from crystal_analysis import Crystal 
from pathlib import Path




cif_file_path = Path(__file__).parent / "../data/test_files/ICSD_CaTiO3_CIF.cif"

cif_file_path = cif_file_path.resolve()

if not cif_file_path.is_file():
    raise FileNotFoundError(f"The CIF file does not exist: {cif_file_path}")

structure = CifParser(cif_file_path).parse_structures(primitive=True)[0]
crystal = Crystal(pymatgen_structure=structure)

eb = crystal.bond_dissociation_enthalpies
vr = crystal.reduction_potentials

print("Type of eb:", type(eb))
print("Contents of eb:", eb)

print("Type of vr:", type(vr))
print("Contents of vr:", vr)