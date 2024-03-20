import pandas as pd
import numpy as np
from pymatgen.core import Composition, Structure
import tarfile
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.bond_valence import BVAnalyzer as BVA
from analysis.crystal_analysis import Crystal
from tqdm import tqdm
from glob import glob

def load_and_merge_csv(data_path: str, csv_pattern: str, merge_columns: list = None, unique_identifier: str = None) -> pd.DataFrame:
    """
    Load CSV files matching a pattern and optionally merge them on specified columns.
    
    Parameters:
    - data_path: Base directory where the CSV files are located.
    - csv_pattern: Pattern to match the CSV files.
    - merge_columns: Columns to merge on if merging multiple CSVs.
    - unique_identifier: Column to add that uniquely identifies rows based on the file they came from.
    
    Returns:
    - A pandas DataFrame containing the merged data.
    """
    csv_paths = sorted(glob(data_path + csv_pattern))
    df_list = [pd.read_csv(csv_path) for csv_path in csv_paths]
    
    if unique_identifier:
        for i, df in enumerate(df_list):
            df[unique_identifier] = i
    df = pd.concat(df_list, ignore_index=True)
    
    if merge_columns:
        df = df.drop_duplicates(subset=merge_columns).reset_index(drop=True)
    
    return df

def filter_defects(df: pd.DataFrame, defect_name: str = "V_O") -> pd.DataFrame:
    """
    Filters the DataFrame for specific defects.
    
    Parameters:
    - df: DataFrame to filter.
    - defect_name: Name of the defect to filter by.
    
    Returns:
    - Filtered DataFrame.
    """
    return df[df["defectname"] == defect_name].reset_index(drop=True)

def add_oxidation_states(df: pd.DataFrame, data_path: str, oxstate_path_pattern: str) -> pd.DataFrame:
    """
    Adds oxidation states to structures within the DataFrame.
    
    Parameters:
    - df: DataFrame containing the structures.
    - data_path: Base directory where the oxidation state files are located.
    - oxstate_path_pattern: Pattern to match the oxidation state files, with placeholders for filename.
    
    Returns:
    - DataFrame with oxidation states added to the structures.

    TODO: consider parallelization?
    """
    for idx, row in df.iterrows():
        oxstate_file = data_path + oxstate_path_pattern.format(filename=row['filename'])
        oxstate = []
        with open(oxstate_file, 'r') as file:
            for x in file.read().split():
                count, state = x.split("*")
                oxstate += int(count) * [float(state)]
        structure = Structure.from_str(row['poscar'], fmt="poscar")
        structure.add_oxidation_state_by_site(oxstate)
        df.at[idx, 'structure'] = structure
    
    return df

def prepare_formula_based_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepares formula-based features like checking if a compound is binary or ternary.
    
    Parameters:
    - df: DataFrame to process.
    
    Returns:
    - DataFrame with additional features.
    """
    df["is_binary_or_ternary"] = df["formula"].apply(lambda x: 2 <= len(Composition(x).elements) <= 3)
    return df

# def calculate_bond_dissociation_energies(df: pd.DataFrame, data_path: str, formula_column: str = "formula") -> pd.DataFrame:
#     """
#     Calculates the sum of bond dissociation energies (Eb) weighted by coordination numbers (CN)
#     for each compound in the dataframe using crystal structure data.

#     Parameters:
#     - df: DataFrame containing the compounds information and paths to their respective CIF files.
#     - data_path: Path to the directory containing CIF files and oxidation state information.
#     - formula_column: Name of the column in df that contains the formula of compounds.

#     Returns:
#     - DataFrame with an additional 'Eb_sum' column containing the calculated sum of bond dissociation energies for each compound.
#     """
#     Eb_sum_list = []
    
#     for _, row in df.iterrows():
#         formula = row[formula_column]
#         cif_path = f"{data_path}/{formula}/structure.cif"  # Modify path as necessary

#         try:
#             structure = Structure.from_file(cif_path)
#             BVA().get_oxi_state_decorated_structure(structure)
#             crystal = Crystal(pymatgen_structure=structure, nn_finder=CrystalNN(weighted_cn=True, cation_anion=True), use_weights=True)
            
#             CN = crystal.cn_dicts
#             Eb = crystal.bond_dissociation_enthalpies
            
#             Eb_sum = sum([np.sum(np.array(list(CN_dict.values())) * np.array(list(Eb_dict.values()))) for CN_dict, Eb_dict in zip(CN, Eb)])
#             Eb_sum_list.append(Eb_sum)
#         except Exception as e:
#             print(f"Error processing {formula}: {e}")
#             Eb_sum_list.append(np.nan)

#     df['Eb_sum'] = Eb_sum_list
#     return df.dropna(subset=['Eb_sum'])


# def calculate_reduction_potentials(df: pd.DataFrame, vr_data_path: str, formula_column: str = "formula") -> pd.DataFrame:
#     """
#     Calculates maximum reduction potentials (Vr) for each compound in the dataframe.

#     Parameters:
#     - df: DataFrame containing the compounds information.
#     - vr_data_path: Path to the CSV file containing Vr values for elements.
#     - formula_column: Name of the column in df that contains the formula of compounds.

#     Returns:
#     - DataFrame with an additional 'vr_max' column containing the maximum Vr for each compound.
#     """
#     vr_df = pd.read_csv(vr_data_path)
#     vr_max_list = []

#     for _, row in df.iterrows():
#         formula = row[formula_column]
#         composition = Composition(formula)
#         metal_vrs = []

#         for element in composition.elements:
#             if element.symbol != "O":
#                 oxi_state_guesses = composition.oxi_state_guesses()
#                 for guess in oxi_state_guesses:
#                     element_vr = vr_df[(vr_df['elem'] == element.symbol) & 
#                                        (vr_df['n'] == guess[element.symbol])]['Vr']
#                     if not element_vr.empty:
#                         metal_vrs.append(element_vr.max())

#         vr_max_list.append(max(metal_vrs) if metal_vrs else np.nan)

#     df['vr_max'] = vr_max_list
#     return df

def calculate_crystal_features(df: pd.DataFrame, structure_column: str = "structure") -> pd.DataFrame:
    """
    Calculates crystal features such as CN-weighted Eb sum and maximum Vr for compounds in the dataframe.

    Parameters:
    - df: DataFrame containing the compounds and their structures.
    - structure_column: Name of the column in df that contains pymatgen Structure objects.

    Returns:
    - DataFrame with additional columns for crystal features ('Eb_sum', 'Vr_max', etc.).
    """
    crystal_features = {'Eb_sum': [], 'Vr_max': []}  # Add more features as necessary

    for _, row in tqdm(df.iterrows(), total=df.shape[0]):
        structure = row[structure_column]
        crystal = Crystal(pymatgen_structure=structure, nn_finder=CrystalNN(weighted_cn=True, cation_anion=True), use_weights=True)

        CN = crystal.cn_dicts
        Eb = crystal.bond_dissociation_enthalpies
        Vr = crystal.reduction_potentials

        # Calculate CN-weighted Eb sum
        Eb_sum = np.sum([np.sum(np.array(list(CN_dict.values())) * np.array(list(Eb_dict.values()))) for CN_dict, Eb_dict in zip(CN, Eb)])
        # Calculate maximum Vr
        Vr_max = max([max(Vr_dict.values()) for Vr_dict in Vr]) if Vr else np.nan

        crystal_features['Eb_sum'].append(Eb_sum)
        crystal_features['Vr_max'].append(Vr_max)

    for feature, values in crystal_features.items():
        df[feature] = values

    return df