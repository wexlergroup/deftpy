from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.core.structure import Structure
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parent.parent.parent))
from analysis.crystal import Crystal
import numpy as np
import pandas as pd
import math
import os
from datetime import datetime

# date_str = datetime.now().strftime("%y%m%d")    # 20241223
# date_str = datetime.now().strftime("%y%b%d")    # 2024Dec23
date_str = datetime.now().strftime("%y-%b-%d")    # 2024-Dec-23

M = 0
N = 18

st = Structure.from_file(f"/Users/manish/PycharmProjects/Proj_CCTM4STCH/2_DFT/1_data/"
                         f"M{M}_N{N}/2_vacancy/0_bulk_structure/CONTCAR_bulk_M{M}_N{N}")

test_cif = f"./junk_{date_str}.cif"


def get_unique_values_from_dicts(dict_list, key):
    seen_pairs = set()
    values = []
    for d in dict_list:
        if key in d:
            pair = (key, d[key])

            if pair not in seen_pairs:
                seen_pairs.add(pair)
                values.append(d[key])
    return values

def get_species_data():
    print()
    unique_values = {}
    Mn_OS = [2, 3, 4]
    Ce_OS = [3, 4]
    species_data = []
    save_file = f"species_data_{date_str}.csv"

    if os.path.exists(save_file):
        df = pd.read_csv(save_file)
        dff = df.drop_duplicates()
        print(dff)
        return dff

    oxi_state_base = {"Ca": 2, "Ti": 4, "O": -2}
    species_template = {"Ca": "Ca{}+", "Ce": "Ce{}+", "Ti": "Ti{}+", "Mn": "Mn{}+"}

    for ceos, mnos in [(ce, mn) for ce in Ce_OS for mn in Mn_OS]:
        print(f"Ce{ceos}, Mn{mnos}")
        oxi_state_override = {**oxi_state_base, "Ce": ceos, "Mn": mnos}
        # print(oxi_state_override)
        st.add_oxidation_state_by_element(oxi_state_override)
        crystal = Crystal(pymatgen_structure=st, nn_finder_type="voronoinn")
        eb = crystal.bond_dissociation_enthalpies
        vr = crystal.reduction_potentials

        for elem in ["Ca", "Ce", "Ti", "Mn"]:
            specie = species_template[elem].format(oxi_state_override[elem])
            try:
                eb_value = next(val for val in get_unique_values_from_dicts(eb, specie))
                vr_value = next(val for val in get_unique_values_from_dicts(vr, specie))
                species_data.append([specie, eb_value, vr_value])
                print(specie)
            except (StopIteration, IndexError):
                print(f"Warning: Data not found for {specie}")

    df = pd.DataFrame(species_data, columns=["species", "eb", "vr"]).sort_values(by="species").drop_duplicates()
    print(df)
    # df.to_csv(save_file, index=False)
    return df

get_species_data()
# print(df)
# get_all()
