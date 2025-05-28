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

# date_str = datetime.now().strftime("%y%m%d")    # 20241222
# date_str = datetime.now().strftime("%y%b%d")    # 2024Dec22
date_str = datetime.now().strftime("%y-%b-%d")    # 2024-Dec-22

M = 0
N = 18

st = Structure.from_file(f"/Users/manish/PycharmProjects/Proj_CCTM4STCH/2_DFT/1_data/"
                         f"M{M}_N{N}/2_vacancy/0_bulk_structure/CONTCAR_bulk_M{M}_N{N}")

cif_junk = f"./test_{date_str}.cif"
def get_all():
    oxi_state_override = {"Ca": 2, "Ce": 3, "Ti": 4, "Mn": 4, "O": -2,}
    st.add_oxidation_state_by_element(oxi_state_override)

    cif_writer = CifWriter(st)
    cif_writer.write_file(cif_junk)
    # cif_file_path = Path(__file__).parent / f"../manish/{cif_junk}"
    cif_file_path = Path(__file__).parent / f"{cif_junk}"
    print(cif_file_path)
    cif_file_path = cif_file_path.resolve()
    if not cif_file_path.is_file():
        raise FileNotFoundError(f"The CIF file does not exist: {cif_file_path}")

    structure = CifParser(cif_file_path).parse_structures(primitive=True)[0]

    crystal = Crystal(pymatgen_structure=structure, nn_finder_type="voronoinn")

    cn = crystal.cn_dicts
    for cdict in cn:
        cdict.pop("O2-", None)

    eb = crystal.bond_dissociation_enthalpies
    for edict in eb:
        edict.pop("O2-", None)

    vr = crystal.reduction_potentials
    for vdict in vr:
        vdict.pop("O2-", None)


    Ce_ox = oxi_state_override["Ce"]
    Ce_specie = f"Ce{Ce_ox}+"
    Mn_ox = oxi_state_override["Mn"]
    Mn_specie = f"Mn{Mn_ox}+"

    print("Type of eb:", type(eb))


    print("Type of vr:", type(vr), len(vr))
    # print("Contents of vr:", vr)

    features = []
    for index, dict in enumerate(cn):
        nCa = dict.get("Ca2+", 4 - dict.get(Ce_specie, 0))
        nCe = dict.get(Ce_specie, 4 - dict.get("Ca2+", 0))
        nTi = dict.get("Ti4+", 2 - dict.get(Mn_specie, 0))
        nMn = dict.get(Mn_specie, 2 - dict.get("Ti4+", 0))

        env = f"{nCe}-{nMn}"

        eb_d = eb[index]
        vr_d = vr[index]
        sigma_eb = (nCa * eb_d.get("Ca2+", 0) + nCe * eb_d.get(Ce_specie, 0) +
                    nTi * eb_d.get("Ti4+", 0) + nMn * eb_d.get(Mn_specie, 0))

        max_vr = max(vr_d.values())
        max_vr_species = max(vr_d, key=vr_d.get)

        features.append([nCa, nCe, nTi, nMn, env, sigma_eb, max_vr, max_vr_species])

    dfff = pd.DataFrame(features, columns=["nCa", "nCe", "nTi", "nMn", "env", "sigma_eb", "max_vr", "max_vr_ion"])
    ddff = dfff.drop_duplicates().sort_values(by="env", ascending=True)
    ddff.to_csv(f"features_sumEb_and_maxVr.csv", index=False)

    print(ddff)


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
    # save_file = "species_data.csv"
    # cif_junk = "./junk.cif"

    # test on Dec 23, 2024
    save_file = "species_data_2024Dec23.csv"
    # cif_junk = "./junk_2025May20.cif"
    if not os.path.exists(save_file):
        for ceos, mnos in zip(Ce_OS, Mn_OS):
            print(f"Ce{ceos}, Mn{mnos}")
            oxi_state_override = {"Ca": 2, "Ce": ceos, "Ti": 4, "Mn": mnos, "O": -2, }
            st.add_oxidation_state_by_element(oxi_state_override)

            cif_writer = CifWriter(st)
            cif_writer.write_file(cif_junk)
            # cif_file_path = Path(__file__).parent / f"../manish/{cif_junk}"
            cif_file_path = Path(__file__).parent / f"{cif_junk}"
            cif_file_path = cif_file_path.resolve()
            if not cif_file_path.is_file():
                raise FileNotFoundError(f"The CIF file does not exist: {cif_file_path}")

            structure = CifParser(cif_file_path).parse_structures(primitive=True)[0]

            crystal = Crystal(pymatgen_structure=structure, nn_finder_type="voronoinn")

            eb = crystal.bond_dissociation_enthalpies
            vr = crystal.reduction_potentials

            Ca_ox = oxi_state_override["Ca"]
            Ca_specie = f"Ca{Ca_ox}+"
            Ce_ox = oxi_state_override["Ce"]
            Ce_specie = f"Ce{Ce_ox}+"
            Ti_ox = oxi_state_override["Ti"]
            Ti_specie = f"Ti{Ti_ox}+"
            Mn_ox = oxi_state_override["Mn"]
            Mn_specie = f"Mn{Mn_ox}+"

            species_list = [Ca_specie, Ce_specie, Ti_specie, Mn_specie]
            for species in species_list:
                eb_value = get_unique_values_from_dicts(eb, species)[0]
                vr_value = get_unique_values_from_dicts(vr, species)[0]
                species_data.append([species, eb_value, vr_value])

        print(species_data)
        df = pd.DataFrame(species_data, columns=["species", "eb", "vr"])
        dff = df.drop_duplicates()
        print(dff)
        # dff.to_csv(save_file, index=False)

    else:
        df = pd.read_csv(save_file)
        dff = df.drop_duplicates()
        print(dff)
        # dff.to_csv(save_file, index=False)

    return dff

df = get_species_data()
