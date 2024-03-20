import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from adjustText import adjust_text

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

    '''
    # set X and y
    X = df_plot.loc[:, ["Eb_sum", "vr_max", "band_gap", 'bv_sum_Crystal', 'o2p_center_from_vbm']]
    y = df_plot.loc['vacancy_formation_energy']

    #get correlations
    corrmap = df_plot.corr()
    top_corr_feat = corrmap.index
    plt.figure(figsize=(8, 8))

    # plt heatmap
    g = sns.heatmap(df_plot[top_corr_feat].corr(), annot=True, annot_kws={"fontsize": 10}, cmap='hot')
    sns.set(font_scale=2)
    labels = [f'E_b', f'V_r', f'E_g', 'BVS', f'O_(2p)']
    plt.xticks(rotation=90)
    plt.yticks(rotation=45)
    plt.tight_layout()

    plt.savefig('heat_map.png')
    plt.show()'''

    #plt BVS from crystal vs nn
    crystal = df_plot['bv_sum_Crystal']
    nearneighbor = df_plot['bv_sum_nn']
    formulas = df_plot['formula']

    deviation = np.abs(crystal - nearneighbor)
    max_dev = np.max(deviation)
    med_dev = np.mean(deviation)

    plt.figure(figsize=(10, 10))
    plt.scatter(crystal, nearneighbor)
    for x, y, dev, formula in zip(crystal, nearneighbor, deviation, formulas):
        if x == y:
            plt.scatter(x, y, color='red', s=50)  # Change color for matching values and adjust point size
        else:
            plt.scatter(x, y, alpha=0.5)  # Reduce alpha for non-matching values to reduce clutter
            if med_dev + 0.2 <= dev <= max_dev:  # Display labels for selected points to reduce clutter
                plt.annotate(formula, (x, y), textcoords="offset points", xytext=(7, 7), ha='left', fontsize=8)
    plt.xlim(-3, 0)
    plt.ylim(-3, 0)

    plt.xlabel('crystal')
    plt.ylabel('near neighbor')
    plt.title('BVS of Crystal vs. Near Neighbor schemes')
    plt.tight_layout()
    plt.savefig('BVS_crystal_vs_nn_colorized.png')
    plt.show()

if __name__ == "__main__":
    main()