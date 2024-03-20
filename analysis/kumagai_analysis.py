import os
import tarfile
import pandas as pd
from tqdm import tqdm
from analysis.crystal_analysis import Crystal
from utils.data_processing import load_and_merge_csv, filter_defects, prepare_formula_based_features, calculate_crystal_features

def main():
    data_path = "../data/papers/kumagai/"

    # load and merge CSV files
    df = load_and_merge_csv(data_path, "oxygen_vacancies_db_data/Li2O/*0.tar.gz", merge_columns=["full_name", "vacancy_index"], unique_identifier="unique_id")

    # remove compounds with transition metals
    df = filter_defects(df, defect_name="V_O")

    # remove unnecessary columns
    df_plot = filter_defects(df_plot, defect_name="V_O")[["formula", "full_name", "band_gap", "formation_energy", "nn_ave_eleneg", "o2p_center_from_vbm", "vacancy_formation_energy", "charge", "vacancy_index"]].reset_index(drop=True)

    # merge bond valence analysis into the df
    df_bva = load_and_merge_csv(data_path, "valence_data_full.csv", merge_columns=["full_name", "vacancy_index"], unique_identifier="unique_id")
    df_plot = pd.merge(df_plot, df_bva[["full_name", "vacancy_index", "valence", "bv_sum_Crystal", "bv_sum_nn"]], on=["full_name", "vacancy_index"], how='outer').drop_duplicates().reset_index(drop=True)

    # prepare formula based features
    df_plot = prepare_formula_based_features(df_plot)

    # calculate crystal features
    df_plot = calculate_crystal_features(df_plot, structure_column="structure")

    # TODO: continue from modeling.py once finished