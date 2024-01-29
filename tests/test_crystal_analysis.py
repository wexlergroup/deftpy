import pytest
from pathlib import Path
from ..analysis.crystal_analysis import Crystal

test_file_dir = Path(__file__).parent.parent / 'data/test_files'

cif_files = list(test_file_dir.glob('*.cif'))
poscar_files = list(test_file_dir.glob('*.txt'))

@pytest.mark.parametrize("cif_file", cif_files)
def test_init(cif_file):
    crystal = Crystal(filepath=str(cif_file))

    assert crystal.structure is not None, f"Initialization failed for {cif_file}"

@pytest.mark.parametrize("poscar_file", poscar_files)
def test_init(poscar_file):
    crystal = Crystal(filepath=str(poscar_file))

    assert crystal.structure is not None, f"Initialization failed for {poscar_file}"
    # all 3 of the CIF files passed, but none of the three OQMD files passed (Element has no attribute oxi_state!)

@pytest.mark.parametrize("input_string, expected_output", [
    ("O2", ["O", "2"]),                 # Element and Quantity
    ("Fe3+", ["Fe", "3+"]),             # Element and Oxidation State
    ("H", ["H"]),                       # Element Symbol Only
    ("Na2O", ["Na", "2O"]),             # Multiple Elements
    ("C60", ["C", "60"]),               # Special Cases (Buckminsterfullerene?)
])
def test_split_before_first_number(input_string, expected_output):
    result = Crystal._split_before_first_number(input_string)

    assert result == expected_output, f"Failed for input '{input_string}'"

def test_parse_species_string():
    pass

def test_initialize_structure_analysis():
    pass

def test_get_values():
    pass

def test_visualize():
    pass


# notes: 