import tarfile
import pandas as pd
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifParser
from glob import glob
import os

#path for scratch slightly different than for analysis scripts
data_path = "/Users/isakov/PycharmProjects/deft_ethans_branch/deftpy/data/papers/kumagai/"

#def NN_from_index(index, structure):

sample_df = pd.DataFrame({
    'formula': ['Al2Ge2O7', 'Al2Ge2O7', 'Ga2O3'],
    'full_name': ['Al2Ge2O7_Va_O1', 'Al2Ge2O7_Va_O2', 'Ga2O3_Va_O3'],
    'band_gap': [3.5, 3.0, 2.2],
    'vacancy_formation_energy': [0.8, 1.2, 1.0],
    'formation_energy': [-3.5, -4.0, -5.2],
    'charge': [2, 0, 3]})

vacancy_indexes = []
for defect in (sample_df["vacancy_formation_energy"].unique()):
    df_defect = sample_df[sample_df["vacancy_formation_energy"] == defect]
    full_name = df_defect["full_name"].iloc[0]
    formula = full_name.split("_")[0]
    vacancy = full_name.split("Va_")[1]
    info_path = os.path.join(data_path, "site_info", formula, "cell_info.txt")
    with open(info_path, "r") as txt:
        index = "Irreducible element: " + vacancy
        lines = txt.readlines()
        for i, line in enumerate(lines):
            if index in line:
                target_line_index = min(i + 5, len(lines) - 1)
                target_line = lines[target_line_index]
                vacancy_index = int(target_line.split("..")[-1])
                vacancy_index = vacancy_index + 1
                print(vacancy_index)
                vacancy_indexes.append(vacancy_index)
                break
print(vacancy_indexes)

cif_path = os.path.join(data_path, "site_info", formula, "supercell.cif")
parser = CifParser(cif_path)
structure = parser.parse_structures(primitive=True)