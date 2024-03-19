import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():

    # define data frames
    df = pd.read_csv('kumagai_Vr_Eb_full_frac.csv')
    df_bva = pd.read_csv("valence_data_full.csv")
    df_bva['vacancy_index'] = df_bva['site'] + 1
    df_bva['merge_on'] = df_bva['full_name'] + df_bva['vacancy_index'].astype(str)
    df_bva = df_bva[["merge_on", 'valence', 'bv_sum_Crystal', 'bv_sum_nn']].reset_index(drop=True)
    df['merge_on'] = df['full_name'] + df['vacancy_index'].astype(str)
    df_2 = pd.merge(df, df_bva, on='merge_on', how='inner')
    df_plot = df_2.drop_duplicates().reset_index(drop=True)
    print(df_plot.columns)
    exit(22)

    # set X and y
    X = df_plot.loc[:, ["Eb_sum", "vr_max", "band_gap", 'bv_sum_Crystal', 'o2p_center_from_vbm']]
    y = df_plot.loc['vacancy_formation_energy']

    #get correlations
    corrmap = df_plot.corr()

if __name__ == "__main__" :
    main()