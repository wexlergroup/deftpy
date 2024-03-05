import pandas as pd
import numpy as np
from pymatgen.analysis.local_env import VoronoiNN, NearNeighbors, CrystalNN
from pymatgen.core import Structure, Composition, IStructure
from pymatgen.analysis.bond_valence import BVAnalyzer as BVA, calculate_bv_sum
import os
from collections import namedtuple

def get_nearest_neighbors(structure, site_index):
    NearestNeighbor = namedtuple('Neighbor', ['site', 'distance'])
    nn = CrystalNN()
    neighbors = nn.get_nn_info(structure, site_index)
    nearest_neighbors = []
    for neighbor_info in neighbors:
        distance = float(neighbor_info['weight'])
        nn_site = neighbor_info['site_index']
        print(nn_site, type(nn_site))
        site = structure[nn_site]
        print(site, type(site))
        nearest_neighbors.append(NearestNeighbor(site, distance))
    return nearest_neighbors

def main():

    df = pd.read_csv(
        '../../../../Library/Application Support/JetBrains/PyCharm2022.2/scratches/complete_df_indexed.csv')
    df["has_transition_metal"] = df.formula.apply(
        lambda x: any([el.is_transition_metal for el in Composition(x)]))
    df = df.loc[~df.has_transition_metal].reset_index(drop=True)
    
    df_plot = df[["formula", "full_name", "band_gap", "formation_energy", "nn_ave_eleneg", "o2p_center_from_vbm",
                      "vacancy_formation_energy", "charge", "vacancy_index"]].reset_index(drop=True)
    valence_for_index = []
    for defect in (df_plot["vacancy_formation_energy"].unique()):
        df_defect = df_plot[df_plot["vacancy_formation_energy"] == defect]
        formula = df_defect['formula'].iloc[0]
        path = '/deftpy/data/papers/kumagai/site_info/'
        cif_path = os.path.join(path, formula, "supercell.cif")
        base_structure = Structure.from_file(cif_path)
        nn_structure = IStructure.from_file(cif_path)

        structure = BVA().get_oxi_state_decorated_structure(base_structure)

        try:
            valences = BVA().get_valences(structure)
            # valences = BVA().get_valences(base_structure) # w/o assigning oxi-states

            # print(valences)
            n = df_defect['vacancy_index'].iloc[0] - 1
            site = structure[n]
            # site = base_structure[n] # w/o assigning oxi-states

            nearest_neighbors = get_nearest_neighbors(structure, n)
            # nearest_neighbors = get_nearest_neighbors(base_structure, n)
            print(nearest_neighbors, type(nearest_neighbors))
            # nn_neighbors = nn_structure.get_neighbors(site=n, r=1)
            # print(type(nearest_neighbors))
            # exit(12)
            # print(type(nn_tuple))
            bv_sum = calculate_bv_sum(site=n, nn_list=nearest_neighbors)
            print(bv_sum)
            exit(12)
            el = site.species_string
            site_valence = valences[n]
            # print(site_valence)
            # print(site, el, site_valence)
            valence_for_index.append({
                "site": n,
                "valence": site_valence,
                "element": el,
            })
        except ValueError:
            valence_for_index.append({
                "site": n,
                "valence": np.nan,
                "element": el,
            })
            pass
    print(valence_for_index)
    
if __name__ == "__main__":
    main()