import os
import tarfile
from glob import glob
import pandas as pd
from pymatgen.analysis.local_env import CrystalNN
from crystal_analysis import Crystal
from pymatgen.core import Structure, Composition
import numpy as np
from tqdm import tqdm
from pymatgen.analysis.bond_valence import BVAnalyzer as BVA

def main():
    data_path = "../data/papers/kumagai/"
    # Li2O test
    with tarfile.open(glob(data_path +"oxygen_vacancies_db_data/Li2O/*0.tar.gz")[0], "r:gz") as tar:
        # Read the member with the name "CONTCAR-finish"
        f = tar.extractfile("CONTCAR-finish")


    df = pd.read_csv("complete_df_indexed.csv")  # aquired using getting_indexes.py

    # *** Remove compounds with transition metals for full set
    df["has_transition_metal"] = df.formula.apply(
        lambda x: any([el.is_transition_metal for el in Composition(x)]))
    df = df.loc[~df.has_transition_metal].reset_index(drop=True)

    # **** Remove unnecessary columns for full data set
    df_plot = df[["formula", "full_name", "band_gap", "formation_energy", "nn_ave_eleneg", "o2p_center_from_vbm",
                  "vacancy_formation_energy", "charge", "vacancy_index"]].reset_index(drop=True)
    # df_plot.to_csv("Kumagai_full.csv")

    # merge bond valence analysis onto the df
    df_bva = pd.read_csv("valence_data_full.csv")
    df_bva['vacancy_index'] = df_bva['site'] + 1
    df_bva['merge_on'] = df_bva['full_name'] + df_bva['vacancy_index'].astype(str)
    df_bva = df_bva[["merge_on", 'valence', 'bv_sum_Crystal', 'bv_sum_nn']].reset_index(drop=True)
    df_plot['merge_on'] = df_plot['full_name'] + df_plot['vacancy_index'].astype(str)
    # print(df_bva['merge_on'], df_plot['merge_on'])
    df_plot2 = pd.merge(df_plot, df_bva, on='merge_on', how='outer')
    df_plot = df_plot2.drop_duplicates()
    # df_plot.to_csv('test.csv')

    # Calculate crystal reduction potentials for binaries and non-binaries - from Vr CSV files
    vr_list = []
    for _, row in df_plot.iterrows():
        # formula = row['formula']
        formula = row["formula"]
        composition = Composition(formula)
        metal_info = []  # For storing information about each metal in the formula
        for el in composition.elements:
            if el.symbol != "O":
                metal = el.symbol
                n_metal = composition.get_el_amt_dict()[metal]
                try:
                    composition = Composition.add_charges_from_oxi_state_guesses(composition)  # guess uses pymatgen get_oxi_state function
                    oxidation_states = Composition.oxi_state_guesses(composition)
                except ValueError:
                    metal_info.append({
                        "metal": metal,
                        "n_metal": n_metal,
                        "oxi_state": np.nan,
                        "vr": np.nan
                    })
                    pass
                try:
                    get_dict = oxidation_states[0]
                    oxi_state = get_dict[metal]
                    # get vr csv as df
                    vr_df = pd.read_csv("../data/features/Vr.csv")
                    # match to metal
                    vr_m_df = vr_df[vr_df["elem"] == metal]
                    # match to oxi_state
                    vr_m_o_df = vr_m_df[vr_m_df["n"] == oxi_state]
                    # choose optimal vr for oxidation
                    vr_max_df = vr_m_o_df[vr_m_o_df['m'] == vr_m_o_df['m'].max()]
                    # print(vr_max_df)
                    vr = vr_max_df["Vr"].iloc[0]
                    # print(vr)
                    # vr = n_atoms * formation_energy / n_metal / oxi_state #original hand calculation of Vr for binary
                    metal_info.append({
                        "metal": metal,
                        "n_metal": n_metal,
                        "oxi_state": oxi_state,
                        "vr": vr
                    })
                except IndexError or ValueError:
                    metal_info.append({
                        "metal": metal,
                        "n_metal": n_metal,
                        "oxi_state": np.nan,
                        "vr": np.nan
                    })
                    pass
        max_vr = max(d['vr'] for d in metal_info)
        vr_list.append(max_vr)

    # Calculate 'vr_max' for the entire DataFrame
    df_plot["vr_max"] = vr_list

    # Save to CSV
    # df_plot.to_csv("kumagai_ternary_vr_from_csv.csv")

    # remove any values where oxi_state could not be assigned
    df_plot = df_plot.dropna()
    # print(df_plot)

    # using corrected structure files
    Eb_sum = []
    for defect in tqdm(df_plot["vacancy_formation_energy"].unique()):
        df_defect = df_plot[df_plot["vacancy_formation_energy"] == defect]
        formula = df_defect["formula"].iloc[0]
        index = df_defect["vacancy_index"].iloc[0]

        with tarfile.open(glob(data_path + "site_info.tar.gz")[0], "r:gz") as tar:
            tar.extractall(path=data_path)
            cif_path = os.path.join(data_path, "site_info", formula, "supercell.cif")
            if os.path.exists(cif_path):
                # MUST use .from_file to ensure the structure is not minimized/primitive or the index will be lost
                base_structure = Structure.from_file(cif_path)
                if base_structure:
                    # for generating oxidation states using Bond Valence (better way - avi)
                    try:
                        structure = BVA().get_oxi_state_decorated_structure(base_structure)
                        
                        crystal = Crystal(pymatgen_structure=structure, n=index, nn_finder=CrystalNN(weighted_cn=True, cation_anion=True), use_weights=True)

                        CN = crystal.cn_dicts
                        Eb = crystal.bond_dissociation_enthalpies


                        for CN_dict, Eb_dict in zip(CN, Eb):
                            CN_array = np.array(list(CN_dict.values()))
                            Eb_array = np.array(list(Eb_dict.values()))
                            Eb_sum.append(np.sum(CN_array * Eb_array))

                    except ValueError:
                        Eb_sum.append(np.nan)
                        print("Eb is nan")
                        pass
                else:
                    print(f"Error: List empty - {base_structure}")
            else:
                print(f"Error: File not found - {cif_path}")
    df_plot["Eb_sum"] = Eb_sum
    df_plot = df_plot.dropna()
    # this data frame will have Vr and Eb_sum in addition to the original columns kept for analysis
    print('complete')
    df_plot.to_csv("kumagai_full_Eb_BV_Vr_from_csv_BVA.csv")

if __name__ == "__main__":
    main()
