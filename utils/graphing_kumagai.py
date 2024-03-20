import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.core import Composition
from sklearn import linear_model
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_score, train_test_split


# Histogram
'''
df_hist = pd.read_csv("kumagai_Vr_Eb_full_frac.csv")
df_hist["element_num"] = df_hist.formula.apply(lambda x:  int(len(Composition(x))))
# plt.hist(df_hist["element_num"], bins=[2, 3, 4, 5], align='mid')
# Counting the frequency of each number of elements
element_count = df_hist['element_num'].value_counts().sort_index()
# Plotting the counts
plt.bar(element_count.index, element_count.values, color='skyblue')
plt.title("distribution of the number of elements in the metal oxide compositions of the Kumagai data set", size=8)
plt.xlabel("number of elements")
plt.xticks([2, 3, 4, 5])
plt.ylabel("number of metal oxides")
plt.savefig("kumagai_full_frac_bar.png")
plt.show()'''
# exit(23)

# # CFM
df_plot = pd.read_csv("kumagai_Vr_Eb_full_frac.csv")
df_bva = pd.read_csv("valence_data_full.csv")
df_bva['vacancy_index'] = df_bva['site'] + 1
df_bva['merge_on'] = df_bva['full_name'] + df_bva['vacancy_index'].astype(str)
df_bva = df_bva[["merge_on", 'valence', 'bv_sum_Crystal', 'bv_sum_nn']].reset_index(drop=True)
df_plot['merge_on'] = df_plot['full_name'] + df_plot['vacancy_index'].astype(str)
df_plot2 = pd.merge(df_plot, df_bva, on='merge_on', how='inner')
df_plot = df_plot2.drop_duplicates().reset_index(drop=True)


# # get binaries from binaries/ternaries
df_plot["is_binary"] = df_plot.formula.apply(lambda x: len(Composition(x)) == 2)

fig, axs = plt.subplots(ncols=3, figsize=(12, 4))

for i, charge in enumerate([0, 1, 2]):
    cfm = linear_model.HuberRegressor()
    # X = df_plot.loc[df_plot.charge == charge, ["vr_max", "band_gap"]]
    # y = df_plot.loc[df_plot.charge == charge, "vacancy_formation_energy"]
    # X = df_plot.loc[df_plot.charge == charge, ["Eb_sum", "vr_max", "band_gap"]]
    # X = df_plot.loc[df_plot.charge == charge, ["Eb_sum", "vr_max", "band_gap", "o2p_center_from_vbm"]]
    # X = df_plot.loc[df_plot.charge == charge, ["Eb_sum", "vr_max", "band_gap", "o2p_center_from_vbm", 'bv_sum_Crystal']]
    X = df_plot.loc[df_plot.charge == charge, ["Eb_sum", "vr_max", "band_gap", 'bv_sum_Crystal']]

    y = df_plot.loc[df_plot.charge == charge, "vacancy_formation_energy"]
    # cfm.fit(X, y)
    cfm.fit(X, y)
    y_pred = cfm.predict(X)
    cfm.score(X, y)

    # model tests
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=21)
    # cfm.fit(X_train, y_train)
    # y_pred = cfm.predict(X_test)
    # mae = mean_absolute_error(y_test, y_pred)

    kf = KFold(n_splits=5, shuffle=True, random_state=21)
    scores = cross_val_score(cfm, X, y, scoring='neg_mean_absolute_error', cv=kf)
    for score in scores:
        print('score for this fold is ', score)
        print('coefficients', cfm.coef_)
    # kf = StratifiedKFold(n_splits=12, shuffle=True, random_state=21)
    # scores = cross_val_score(cfm, X, y, scoring='neg_mean_absolute_error', cv=12)
    # scores = cross_val_score(cfm, X_train, y_train, scoring='neg_mean_absolute_error', cv=kf)
    mean = np.mean(scores)
    std = np.std(scores)
    print('for model, mean score = ', mean)
    rmse = mean_squared_error(y, y_pred, squared=False)
    # rmse = mean_squared_error(y_test, y_pred, squared=False)
    confidence = (np.quantile(scores, [0.025, .975]))

    # define colors for binary, ternary etc

    # define colors for element groups
    df_colors = df_plot[df_plot["charge"] == i]
    periodic_table_groups = {
        "H": 1, "He": 18,
        "Li": 1, "Be": 2, "B": 13, "C": 14, "N": 15, "O":16, "F": 17, "Ne": 18,
        "Na": 1, "Mg": 2, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
        "K": 1, "Ca": 2, "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12,
        "Ga": 13, "Ge": 14, "As": 15, "Se": 16, "Br": 17, "Kr": 18,
        "Rb": 1, "Sr": 2, "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10, "Ag": 11, "Cd": 12,
        "In": 13, "Sn": 14, "Sb": 15, "Te": 16, "I": 17, "Xe": 18,
        "Cs": 1, "Ba": 2, "La": 3, "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12,
        "Tl": 13, "Pb": 14, "Bi": 15, "Po": 16, "At": 17, "Rn": 18
    }
    # color_map = {1: 'red', 2: 'blue', 3: 'green', 4: 'orange', 5: 'purple', 6: 'brown', 7: 'pink', 8: 'gray', 9: 'olive', 10: 'cyan', 11: 'magenta', 12: 'yellow', 13: 'lime', 14: 'teal', 15: 'navy', 16: 'maroon', 17: 'skyblue', 18: 'black'}
    # color map w/o TM
    color_map = {1: 'red', 2: 'blue', 13: 'lime', 14: 'teal', 15: 'navy', 16: 'maroon', 17: 'skyblue', 18: 'black'}
    shape_map = {1: '*', 2: 's', 13: 'h', 14: 'x', 15: 'D', 16: '.', 17: 'o', 18: 'v'}
    def get_group(element):
        if element in periodic_table_groups:
            return periodic_table_groups[element]
        else:
            return None

    # df_colors['Groups'] = df_colors['formula'].str.findall(r'[A-Z][a-z]*').apply(lambda x: [get_group(e) for e in x])
    # df_colors['group1'] = df_colors['Groups'].apply(lambda x: x[0])
    # df_colors['group2'] = df_colors['Groups'].apply(lambda x: x[1])
    # # print(df_colors['group1'])
    # colors = df_colors['group1'].map(color_map)
    # markers = df_colors['group2'].map(shape_map)
    #
    # for yy, pred, cc, mm in zip(y, y_pred, colors, markers):
    #     axs[i].scatter(pred, yy, marker=mm)
    # plt.show()
    # exit(23)
    # # print(color_map)
    # # exit(32)

    # Plot results
    axs[i].plot(y_pred, y, "o")
    # axs[i].scatter(y_pred, y, c=colors)
    # axs[i].scatter(y_pred, y_test)
    # axs[i].scatter(y_pred, y, c=colors, marker=markers)

    # Plot parity line
    # axs[i].plot([-4, 10], [-4, 10], "--", color="black")
    axs[i].plot([-4, 10], [-4, 10], "--")

    # Set axis limits
    axs[i].set_xlim(-4, 10)
    axs[i].set_ylim(-4, 10)

    # Add equation
    # equation = "$E_v = {:.2f} {:+.2f} E_b {:+.2f} V_r {:+.2f} E_g$".format(cfm.intercept_, cfm.coef_[0], cfm.coef_[1],
    #                                                                        cfm.coef_[2])
    # equation = "$E_v = {:.2f} {:+.2f} E_b {:+.2f} V_r {:+.2f} E_g {:+.2f} O_p$".format(cfm.intercept_, cfm.coef_[0], cfm.coef_[1], cfm.coef_[2], cfm.coef_[3])
    # equation = "$E_v = {:.2f} {:+.2f} E_b {:+.2f} V_r {:+.2f} E_g {:+.2f} O_2p {:+.2f} BV_sum$".format(cfm.intercept_, cfm.coef_[0], cfm.coef_[1], cfm.coef_[2], cfm.coef_[3], cfm.coef_[4])
    equation = "$E_v = {:.2f} {:+.2f} E_b {:+.2f} V_r {:+.2f} E_g {:+.2f} BVS$".format(cfm.intercept_, cfm.coef_[0], cfm.coef_[1], cfm.coef_[2], cfm.coef_[3])


    axs[i].set_xlabel(equation)

    # Add MAE
    mae = mean_absolute_error(y, y_pred)
    # mae = mean_absolute_error(y_test, y_pred)
    axs[i].text(0.1, 0.9, "MAE = {:.2f} eV".format(mae), size=9, transform=axs[i].transAxes)

    # Add number of data points
    axs[i].text(0.1, 0.85, f"n = {len(y)}", size=9, transform=axs[i].transAxes)

    #ADD KFold score analysis
    axs[i].text(0.1, 0.8, "KFold = 5", size=9, transform=axs[i].transAxes)
    axs[i].text(0.1, 0.75, f"mean score = {mean:.2f}", size=9, transform=axs[i].transAxes)
    axs[i].text(0.1, 0.70, f"stand. dev. = {std:.2f}", size=9, transform=axs[i].transAxes)
    axs[i].text(0.1, 0.65, f"RMSE = {rmse:.2f}", size=9, transform=axs[i].transAxes)
    # axs[i].text(0.1, 0.55, f"95% interval = {confidence:.2f}", size=9, transform=axs[i].transAxes)
    # Add charge as title
    axs[i].set_title(f"Charge {charge}")

    # Add Legend
    # legends = {"binary": "blue", "ternary": "red"}
    # axs[i].legend(colors, loc="lower right")

    # Add y-axis label
    if i == 0:
        axs[i].set_ylabel("$E_v$ (eV)")
    # if i == 2:
    #     for x in color_map:
    #         axs[i].scatter([], [], c=color_map[x], label=f'group {x}')
    #     axs[i].legend(loc='lower right', bbox_to_anchor=(1.5, 0.15))
        # for x in shape_map:
        #     axs[i].scatter([], [], marker=shape_map[x], label=f'group {x}')
        # axs[i].legend(loc='lower left', bbox_to_anchor=(2, 0.15))
plt.tight_layout()
# plt.savefig("kumagai_full_vr_eb_frac_K5_BVA.png", dpi=300)
plt.show()
