import csv

import pandas as pd
import numpy as np
from pymatgen.analysis.local_env import VoronoiNN, NearNeighbors, CrystalNN
from pymatgen.core import Structure, Composition, IStructure, Element
from pymatgen.analysis.bond_valence import BVAnalyzer as BVA, calculate_bv_sum
import os
from tqdm import tqdm
from collections import namedtuple

def get_nearest_neighbors(structure, site_index):
    # NearestNeighbor = namedtuple('Neighbor', ['Site'])
    # nn = CrystalNN()
    nn = CrystalNN()
    neighbors = nn.get_nn_info(structure, site_index)
    nearest_neighbors = []
    for neighbor_info in neighbors:
        distance = float(neighbor_info['weight'])
        nn_site_index = neighbor_info['site_index']
        site = neighbor_info['site']
        # site = structure[nn_site]
        nearest_neighbors.append(site)
        # nn_site = neighbor_info['site']
        # print(type(nn_site))
        # species = structure[nn_site_index]
        # print(type(species))
        # nearest_neighbors.append((NearestNeighbor(nn_site, distance, species)))
    return nearest_neighbors

def main():

    df = pd.read_csv('complete_df_indexed.csv')
    df["has_transition_metal"] = df.formula.apply(lambda x: any([el.is_transition_metal for el in Composition(x)]))
    df = df.loc[~df.has_transition_metal].reset_index(drop=True)
    
    df_plot = df[["formula", "full_name", "band_gap", "formation_energy", "nn_ave_eleneg", "o2p_center_from_vbm",
                      "vacancy_formation_energy", "charge", "vacancy_index"]].reset_index(drop=True)
    valence_for_index = []
    for defect in tqdm(df_plot["vacancy_formation_energy"].unique()):
        df_defect = df_plot[df_plot["vacancy_formation_energy"] == defect]
        formula = df_defect['formula'].iloc[0]
        path = '../data/papers/kumagai/site_info/'
        cif_path = os.path.join(path, formula, "supercell.cif")
        base_structure = Structure.from_file(cif_path)
        nn_structure = IStructure.from_file(cif_path)

        try:
            # valences = BVA().get_valences(structure)
            valences = BVA().get_valences(base_structure) # w/o assigning oxi-states

            # print(valences)
            n = df_defect['vacancy_index'].iloc[0] - 1
            # defect_site = structure[n]
            defect_site = base_structure[n] # w/o assigning oxi-states
            try:
                structure = BVA().get_oxi_state_decorated_structure(base_structure)
                nearest_neighbors = get_nearest_neighbors(structure, n)
            except ValueError:
                nearest_neighbors = get_nearest_neighbors(base_structure, n)
            pass
            # print(nearest_neighbors, type(nearest_neighbors))
            Neighbors = base_structure.get_neighbors(base_structure[n], 3)
            # print(Neighbors, type(Neighbors))
            bv_sum_defined = calculate_bv_sum(site=defect_site, nn_list=nearest_neighbors)
            bv_sum_Neighbors = calculate_bv_sum(site=defect_site, nn_list=Neighbors)
            # print(bv_sum_Neighbors, bv_sum_defined)
            el = defect_site.species_string
            site_valence = valences[n]

            valence_for_index.append({
                "site": n,
                "formula": formula,
                "valence": site_valence,
                "bv_sum_Crystal": bv_sum_defined,
                "bv_sum_nn": bv_sum_Neighbors,
                "element": el,
            })
        except ValueError:
            valence_for_index.append({
                "site": n,
                "formula": formula,
                "valence": np.nan,
                "bv_sum_Crystal": bv_sum_defined,
                "bv_sum_nn": bv_sum_Neighbors,
                "element": el,
            })
            pass
    csv_file_path = "valence_data.csv"

    # Define the field names
    field_names = ["site", "formula", "valence", "bv_sum_Crystal", "bv_sum_nn", "element"]

    # Write the data to the CSV file
    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=field_names)

        # Write the header
        writer.writeheader()

        # Write the data
        for data_row in valence_for_index:
            writer.writerow(data_row)
    # print(valence_for_index)
    
if __name__ == "__main__":
    main()