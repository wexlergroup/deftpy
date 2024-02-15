""" Get crystal features for structures in Yu Kumagai's Physical Review Materials Paper """
import os
import tarfile
import warnings
from glob import glob

import matplotlib.pyplot as plt
import pandas as pd
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from sklearn import linear_model
from sklearn.metrics import mean_absolute_error
from crystal_analysis import Crystal
from pymatgen.core import Structure, Composition, Element
import subprocess
import numpy as np
from tqdm import tqdm
from pymatgen.analysis.bond_valence import BVAnalyzer as BVA
def main():
    data_path = "../data/papers/kumagai/"
    # Li2O
    with tarfile.open(glob(data_path +"oxygen_vacancies_db_data/Li2O/*0.tar.gz")[0], "r:gz") as tar:
        # Read the member with the name "CONTCAR-finish"
        f = tar.extractfile("CONTCAR-finish")
        #print(f.read().decode("utf-8"))

    # Get data
    # df_0 = pd.read_csv(data_path + "vacancy_formation_energy_ml/charge0.csv")  # neutral vacancies
    # df_1 = pd.read_csv(data_path + "vacancy_formation_energy_ml/charge1.csv")  # +1 charged vacancies
    # df_2 = pd.read_csv(data_path + "vacancy_formation_energy_ml/charge2.csv")  # +2 charged vacancies

    # Add charge column
    # df_0["charge"] = 0
    # df_1["charge"] = 1
    # df_2["charge"] = 2

    # Combine dataframes
    # df = pd.concat([df_0, df_1, df_2], ignore_index=True).reset_index(drop=True)

    # Remove the column named "Unnamed: 0"
    # df = df.drop("Unnamed: 0", axis=1)

    df = pd.read_csv("complete_df_indexed.csv")  # aquired using getting_indexes.py

    # Remove non-binary compounds
    # df["is_binary"] = df.formula.apply(lambda x: len(Composition(x)) == 2)
    # df_binary = df.loc[df.is_binary].reset_index(drop=True)
    #
    # # Remove compounds with transition metals
    # df_binary["has_transition_metal"] = df_binary.formula.apply(
    #     lambda x: any([el.is_transition_metal for el in Composition(x)]))
    # df_binary = df_binary.loc[~df_binary.has_transition_metal].reset_index(drop=True)

    # **Remove compounds w/ transition metals including ternary compounds
    df["is_binary_or_ternary"] = df["formula"].apply(lambda x: 2 <= len(Composition(x).elements) <= 3)
    df_ternary = df.loc[df.is_binary_or_ternary].reset_index(drop=True)

    df_ternary["has_transition_metal"] = df_ternary.formula.apply(
        lambda x: any([el.is_transition_metal for el in Composition(x)]))
    df_ternary = df_ternary.loc[~df_ternary.has_transition_metal].reset_index(drop=True)

    # *** Remove compounds with transition metals for full set
    # df["has_transition_metal"] = df.formula.apply(
    #     lambda x: any([el.is_transition_metal for el in Composition(x)]))
    # df = df.loc[~df.has_transition_metal].reset_index(drop=True)

    # Remove unnecessary columns
    # df_plot = df_binary[["formula", "full_name", "band_gap", "formation_energy", "nn_ave_eleneg", "o2p_center_from_vbm",
    #                      "vacancy_formation_energy", "charge"]].reset_index(drop=True)
    # df_plot.to_csv("Kumagai_binary_clean.csv")
    # exit(4)

    # ** Remove unnecessary columns including non-binary compounds
    # To return to the binaries this line must be commented out
    df_plot = df_ternary[["formula", "full_name", "band_gap", "formation_energy", "nn_ave_eleneg", "o2p_center_from_vbm",
                         "vacancy_formation_energy", "charge", "vacancy_index"]].reset_index(drop=True)

    # df_plot.to_csv("Kumagai_ternary.csv")

    # **** Remove unnecessary columns for full data set
    # df_plot = df[["formula", "full_name", "band_gap", "formation_energy", "nn_ave_eleneg", "o2p_center_from_vbm",
    #                     "vacancy_formation_energy", "charge"]].reset_index(drop=True)
    # df_plot.to_csv("Kumagai_full.csv")

    # Calculate crystal reduction potentials for binaries
    '''n_atoms = []  # number of atoms in the compound formula
    metals = []  # metal in the compound formula
    n_metals = []  # number of metal atoms in the compound formula
    oxi_states = []  # oxidation state of the metal in the compound formula
    for _, row in df_plot.iterrows():
        formula = row["formula"]
        composition = Composition(formula)
        metal = [el.symbol for el in composition.elements if el.symbol != "O"][0]
        n_metal = composition.get_el_amt_dict()[metal]
        n_oxygen = composition.get_el_amt_dict()["O"]
        oxi_state = 2 * n_oxygen / n_metal
        n_atoms.append(composition.num_atoms)
        metals.append(metal)
        n_metals.append(n_metal)
        oxi_states.append(oxi_state)
    df_plot["n_atoms"] = n_atoms
    df_plot["metal"] = metals
    df_plot["n_metal"] = n_metals
    df_plot["oxi_state"] = oxi_states
    df_plot["vr"] = df_plot["n_atoms"] * df_plot["formation_energy"] / df_plot["n_metal"] / df_plot["oxi_state"]
    df_plot["vr_max"] = df_plot["vr"].max()
    df_plot.to_csv("Kumagai_binaries.csv")
    #print(df_plot[["formula", "metal", "oxi_state"]])'''

    # Calculate crystal reduction potentials including non-binaries
    vr_list = []
    for _, row in df_plot.iterrows():
        formula = row["formula"]
        formation_energy = row["formation_energy"]
        composition = Composition(formula)
        n_atoms = composition.num_atoms
        metal_info = []  # For storing information about each metal in the formula
        elements_list = [el.symbol for el in composition.elements]
        for el in composition.elements:
            if el.symbol != "O":
                metal = el.symbol
                # print(metal)
                n_metal = composition.get_el_amt_dict()[metal]
                # print(n_metal)
                n_oxygen = composition.get_el_amt_dict()["O"]
                # print(n_oxygen)
                '''if len(composition) != 2:
                    element_amount_dict = composition.get_el_amt_dict()
                    print(element_amount_dict)

                    n_alt = [element_amount_dict[element] for element in element_amount_dict if element != "O" and element != metal][0]
                    if metal in ["Li", "Na", "K", "Rb", "Cs"]:
                        oxi_state = +1.0
                    elif metal in ["Be", "Mg", "Ca", "Sr", "Ba"]:
                        oxi_state = +2.0
                    elif [el for el in elements_list if el != "O" and el != metal] in ["Li", "Na", "K", "Rb", "Cs"]:
                        oxi_alt = +1.0
                        oxi_state = (2 * n_oxygen - oxi_alt * n_alt) / n_metal
                    elif [el for el in elements_list if el != "O" and el != metal] in ["Be", "Mg", "Ca", "Sr", "Ba"]:
                        oxi_alt = +2.0
                        oxi_state = (2 * n_oxygen - oxi_alt * n_alt) / n_metal
                    else:
                        oxi_state = 2 * n_oxygen / (n_metal + n_alt)
                else:
                    oxi_state = 2 * n_oxygen / n_metal'''
                try:
                    composition = Composition.add_charges_from_oxi_state_guesses(composition)
                    # print(composition)
                    oxidation_states = Composition.oxi_state_guesses(composition)
                    # print(oxidation_states)
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
                    # print(oxi_state)
                    # print(metal)
                    # get vr csv as df
                    vr_df = pd.read_csv("../data/features/Vr.csv")
                    # match to metal
                    vr_m_df = vr_df[vr_df["elem"] == metal]
                    # print(vr_m_df)
                    # match to oxi_state
                    vr_m_o_df = vr_m_df[vr_m_df["n"] == oxi_state]
                    # print(vr_m_o_df)
                    # choose optimal vr for oxidation
                    vr_max_df = vr_m_o_df[vr_m_o_df['m'] == vr_m_o_df['m'].max()]
                    # print(vr_max_df)
                    vr = vr_max_df["Vr"].iloc[0]
                    # print(vr)
                    # vr = n_atoms * formation_energy / n_metal / oxi_state
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
        # print(metal_info)
        max_vr = max(d['vr'] for d in metal_info)
        vr_list.append(max_vr)
    # print(len(vr_list))
    # print(len(oxi_states))
    # print(len(df_plot))

    # Calculate 'vr_max' for the entire DataFrame
    df_plot["vr_max"] = vr_list
    # df_plot["oxi_state"] = oxi_states

    # Save to CSV
    # df_plot.to_csv("kumagai_ternary_vr_from_csv.csv")

    # remove any values where oxi_state could not be assigned
    df_plot = df_plot.dropna()

    # using corrected structure files
    structures = []
    Eb_sum = []

    for defect in tqdm(df_plot["vacancy_formation_energy"].unique()):
        df_defect = df_plot[df_plot["vacancy_formation_energy"] == defect]
        # print(df_defect)
        formula = df_defect["formula"].iloc[0]
        # print(formula)
        full_name = df_defect["full_name"].iloc[0]
        # print(full_name)
        # metal = df_defect["metal"].iloc[0]
        # print(metal)
        index = df_defect["vacancy_index"].iloc[0]

        with tarfile.open(glob(data_path + "site_info.tar.gz")[0], "r:gz") as tar:
            tar.extractall(path=data_path)
            cif_path = os.path.join(data_path, "site_info", formula, "supercell.cif")
            # print(cif_path)
            if os.path.exists(cif_path):
                    base_structure = Structure.from_file(cif_path)
                # print(strucutre)
                # exit(3)
                # parser = CifParser(cif_path)
                # parser_list = parser.parse_structures(primitive=True)
                # print(f"parser list length: {len(parser_list)}")
                # if parser_list:
                #     structure = parser.parse_structures(primitive=True)[0]
                    if base_structure:
                        '''oxi_states = {"O": -2, metal: oxi_state_metal}
                        # structure_copy = structure.copy()
                        structure.add_oxidation_state_by_element(oxidation_states=oxi_states)
                        #print("ox states added")
                        structures.append(structure)'''

                        # potential alternative for oxidation states ... better way - avi
                        try:
                            structure = BVA().get_oxi_state_decorated_structure(base_structure)

                            # # if you want to check symmetries/indexing
                            # analyzer = SpacegroupAnalyzer(structure)
                            # symmetry_dataset = analyzer.get_symmetry_dataset()
                            #
                            # wyckoff_symbols = symmetry_dataset['wyckoffs']
                            # wyckoff = wyckoff_symbols[index]
                            #
                            # print(formula, index, f"wyckoff = {wyckoff}")
                            # # structures.append(structure)

                            crystal = Crystal(pymatgen_structure=structure, n=index)
                            # print("structure passed to crystal object")

                            CN = crystal.cn_dicts
                            Eb = crystal.bond_dissociation_enthalpies

                            # Eb_sum = [
                            #     np.sum(np.array(list(cn.values())) * np.array(list(be.values())))
                            #     for cn, be in zip(CN, Eb)
                            # ]

                            # print(f"Eb_sum length: {len(Eb_sum)}")

                            for CN_dict, Eb_dict in zip(CN, Eb):
                                CN_array = np.array(list(CN_dict.values()))
                                Eb_array = np.array(list(Eb_dict.values()))
                                Eb_sum.append(np.sum(CN_array * Eb_array))
                                # print(CN_array, Eb_array)
                            # print("length is " + str(len(Eb_sum)))
                        except ValueError:
                            Eb_sum.append(np.nan)
                            pass
                    else:
                        print(f"Error: List empty - {base_structure}")
                        # print(f"Error: List empty - {parser_list}")
            else:
                print(f"Error: File not found - {cif_path}")
    df_plot["Eb_sum"] = Eb_sum
    df_plot = df_plot.dropna()
    df_plot.to_csv("kumagai_Eb_Vr.csv")
    exit(234)

    #calculate sum Eb
    '''structures = []
    Eb_sum = []
    for defect in tqdm(df_plot["vacancy_formation_energy"].unique()):
    #for full_name in tqdm(df_plot["full_name"].unique()):
        # these four lines might be unnecessary ... idrk
        df_defect = df_plot[df_plot["vacancy_formation_energy"] == defect]
        full_name = df_defect["full_name"].iloc[0]
        print(full_name)
        formula = df_defect["formula"].iloc[0]
        charge = df_plot["charge"].iloc[0]

        # Open the outer tar.gz file
        with tarfile.open(glob(data_path + "oxygen_vacancies_db_data/" + formula + ".tar.gz")[0], "r:gz") as outer_tar:
            # Specify the path to the inner tar.gz file within the outer tar.gz file
            inner_tar_path = str(str(formula) + "/" + str(full_name) + "_" + str(charge) + ".tar.gz")

            # Extract the inner tar.gz file from the outer tar.gz file
            inner_tar_info = outer_tar.getmember(inner_tar_path)
            inner_tar_file = outer_tar.extractfile(inner_tar_info)

            # Open the inner tar.gz file
            with tarfile.open(fileobj=inner_tar_file, mode="r:gz") as inner_tar:
                # obtain the contcar file
                d = inner_tar.extractfile("CONTCAR-finish")
                contcar = d.read().decode("utf-8")
                poscar = Poscar.from_str(contcar)
                #crystal = Crystal(poscar_string=poscar)
                #pass poscar to a structure object
                structure = poscar.structure
                #assign oxidation states
                oxi_states = {"O": -2, str(df_plot.loc[df_plot['full_name'] == full_name, "metal"].iloc[0]): float(df_plot.loc[df_plot['full_name'] == full_name, "oxi_state"].iloc[0])}
                structure_copy = structure.copy()
                structure_copy.add_oxidation_state_by_element(oxidation_states=oxi_states)
                #print(structure_copy)
                crystal = Crystal(pymatgen_structure=structure_copy)
                structures.append(crystal)
                #print(structures)
                CN = crystal.cn_dicts
                Eb = crystal.bond_dissociation_enthalpies
                #Vr = crystal.reduction_potentials

                # Calculate CN-weighted Eb sum

                for CN_dict, Eb_dict in zip(CN, Eb):
                    CN_array = np.array(list(CN_dict.values()))
                    Eb_array = np.array(list(Eb_dict.values()))
                    Eb_sum.append(np.sum(CN_array * Eb_array))

    print(Eb_sum)
    exit(12)'''

    # Fit basic crystal feature model (cfm)
    fig, axs = plt.subplots(ncols=3, figsize=(12, 4))
    for i, charge in enumerate([0, 1, 2]):
        cfm = linear_model.HuberRegressor()
        X = df_plot.loc[df_plot.charge == charge, ["vr_max", "band_gap"]]
        y = df_plot.loc[df_plot.charge == charge, "vacancy_formation_energy"]
        cfm.fit(X, y)
        y_pred = cfm.predict(X)

        # Plot results
        axs[i].plot(y_pred, y, "o")

        # Plot parity line
        axs[i].plot([-4, 10], [-4, 10], "--")

        # Set axis limits
        axs[i].set_xlim(-4, 10)
        axs[i].set_ylim(-4, 10)

        # Add equation
        equation = "$E_v = {:.2f} {:+.2f} V_r {:+.2f} E_g$".format(cfm.intercept_, cfm.coef_[0], cfm.coef_[1])
        axs[i].set_xlabel(equation)

        # Add MAE
        mae = mean_absolute_error(y, y_pred)
        axs[i].text(0.1, 0.9, "MAE = {:.2f} eV".format(mae), size=9, transform=axs[i].transAxes)

        # Add number of data points
        axs[i].text(0.1, 0.8, f"n = {len(y)}", size=9, transform=axs[i].transAxes)

        # Add charge as title
        axs[i].set_title(f"Charge {charge}")

        # Add y-axis label
        if i == 0:
            axs[i].set_ylabel("$E_v$ (eV)")

    plt.tight_layout()
    plt.savefig("kumagai_ternary_vr_from_csv.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()