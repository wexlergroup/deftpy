import pytest
import pandas as pd
from pathlib import Path
from pymatgen.core import Species, Structure
from unittest.mock import Mock, patch
from ..analysis.crystal_analysis import Crystal

test_file_dir = Path(__file__).parent.parent / 'data/test_files'

cif_files = list(test_file_dir.glob('*.cif'))
poscar_files = list(test_file_dir.glob('*.txt'))

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


@pytest.mark.parametrize("input_string, expected_species, expected_symbol, expected_state", [
    ("Fe2+", Species("Fe", 2), "Fe", 2),       # Valid with oxidation state
    ("O", None, "O", 0),                       # Valid without oxidation state
    ("InvalidString", None, "InvalidString", 0),  # Invalid format
    ("", None, "", 0),                         # Empty string
])
def test_parse_species_string(input_string, expected_species, expected_symbol, expected_state):
    species, symbol, state = Crystal._parse_species_string(input_string)

    assert species == expected_species, f"Failed for input '{input_string}'"
    assert symbol == expected_symbol, f"Failed for input '{input_string}'"
    assert state == expected_state, f"Failed for input '{input_string}'"

@pytest.mark.parametrize("cif_file", cif_files)
def test_initialize_cif_analysis(cif_file):
    valid_structure = cif_file
    crystal = Crystal(filepath=str(valid_structure))
    cn_dicts = crystal._initialize_structure_analysis()

    assert cn_dicts is not None, "cn_dicts should be initialized"
    assert isinstance(cn_dicts, list), "cn_dicts should be a list"

    cn_dicts_second_call = crystal._initialize_structure_analysis()
    assert cn_dicts == cn_dicts_second_call, "Repeated calls should yield the same result"

@pytest.mark.parametrize("poscar_file", poscar_files)
def test_initialize_poscar_analysis(poscar_file):
    
    valid_structure = poscar_file
    crystal = Crystal(filepath=str(valid_structure))
    cn_dicts = crystal._initialize_structure_analysis()

    assert cn_dicts is not None, "cn_dicts should be initialized"
    assert isinstance(cn_dicts, list), "cn_dicts should be a list"

    cn_dicts_second_call = crystal._initialize_structure_analysis()
    assert cn_dicts == cn_dicts_second_call, "Repeated calls should yield the same result"

# def test_initialize_structure_analysis_invalid_structure():
#     invalid_structure = None

#     crystal = Crystal(pymatgen_structure=invalid_structure)

#     with pytest.raises(ValueError):
#         crystal._initialize_structure_analysis()

def test_get_values():
    crystal = Crystal(filepath=str(cif_files[0]))

    eb_values = crystal._get_values(crystal.eb, "Eb", "os")
    vr_values = crystal._get_values(crystal.vr, "Vr", "n")

    assert isinstance(eb_values, list), "Returned eb_values should be a list"
    assert isinstance(vr_values, list), "Returned vr_values should be a list"