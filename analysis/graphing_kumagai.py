import pandas as pd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.core import Composition
from pymatgen.core.periodic_table import ElementBase
from sklearn import linear_model
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import KFold, cross_val_score
from mendeleev import element

# Histogram

# df_hist = pd.read_csv("kumagai_Vr_Eb_full_frac.csv")
# df_hist["element_num"] = df_hist.formula.apply(lambda x:  int(len(Composition(x))))
# plt.hist(df_hist["element_num"], align='mid')
# plt.title("distribution of elements in metal oxide compositions")
# plt.xlabel("number of elements")
# plt.ylabel("number of metal oxides")
# plt.savefig("kumagai_full_frac_histogram.png")
# plt.show()

# CFM
df_plot = pd.read_csv("kumagai_Eb_Vr_frac.csv")

# # get binaries from binaries/ternaries
# df_plot["is_binary"] = df_plot.formula.apply(lambda x: len(Composition(x)) == 2)

fig, axs = plt.subplots(ncols=3, figsize=(12, 4))

for i, charge in enumerate([0, 1, 2]):
    cfm = linear_model.HuberRegressor()
    #X = df_plot.loc[df_plot.charge == charge, ["vr_max", "band_gap"]]
    #y = df_plot.loc[df_plot.charge == charge, "vacancy_formation_energy"]
    X = df_plot.loc[df_plot.charge == charge, ["Eb_sum", "vr_max", "band_gap"]]
    y = df_plot.loc[df_plot.charge == charge, "vacancy_formation_energy"]
    cfm.fit(X, y)
    y_pred = cfm.predict(X)

    # model tests
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=21)
    # cfm.fit(X_train, y_train)
    # y_pred = cfm.predict(X_test)
    # mae = mean_absolute_error(y_test, y_pred)

    kf = KFold(n_splits=5, shuffle=True, random_state=21)
    scores = cross_val_score(cfm, X, y, scoring='neg_mean_absolute_error', cv=kf)
    mean = np.mean(scores)
    std = np.std(scores)
    rmse = mean_squared_error(y, y_pred, squared=False)
    confidence = (np.quantile(scores, [0.025, .975]))

    # define colors for binary, ternary etc
    # df_colors = df_plot[df_plot["charge"] == i]
    # colors = {True: 'blue', False: 'red'}
    # color_map = df_colors["is_binary"].map(colors)
    # df_colors = df_plot[df_plot["charge"] == i]
    # df_colors["element_num"] = df_colors.formula.apply(lambda x: int(len(Composition(x))))
    # colors = {2: 'blue', 3: 'red', 4: 'green', 5: 'orange'}
    # # color_map = df_colors["element_num"].map(colors)

    # define colors for element groups
    # df_colors = df_plot[df_plot["charge"] == i]
    # periodic_table_groups = {
    #     "H": 1, "He": 18,
    #     "Li": 1, "Be": 2, "B": 13, "C": 14, "N": 15, "O":16, "F": 17, "Ne": 18,
    #     "Na": 1, "Mg": 2, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    #     "K": 1, "Ca": 2, "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12,
    #     "Ga": 13, "Ge": 14, "As": 15, "Se": 16, "Br": 17, "Kr": 18,
    #     "Rb": 1, "Sr": 2, "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10, "Ag": 11, "Cd": 12,
    #     "In": 13, "Sn": 14, "Sb": 15, "Te": 16, "I": 17, "Xe": 18,
    #     "Cs": 1, "Ba": 2, "La": 3, "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12,
    #     "Tl": 13, "Pb": 14, "Bi": 15, "Po": 16, "At": 17, "Rn": 18
    # }
    # color_map = {1: 'red', 2: 'blue', 3: 'green', 4: 'orange', 5: 'purple', 6: 'brown', 7: 'pink', 8: 'gray', 9: 'olive', 10: 'cyan', 11: 'magenta', 12: 'yellow', 13: 'lime', 14: 'teal', 15: 'navy', 16: 'maroon', 17: 'skyblue', 18: 'black'}
    #
    # def get_group(element):
    #     if element in periodic_table_groups:
    #         return periodic_table_groups[element]
    #     else:
    #         return None
    #
    # df_colors['Groups'] = df_colors['formula'].str.findall(r'[A-Z][a-z]*').apply(lambda x: [get_group(e) for e in x])
    # df_colors['group'] = df_colors['Groups'].apply(lambda x: x[0])
    # # print(df_colors['group'])
    # colors = df_colors['group'].map(color_map)
    # # print(color_map)
    # # exit(32)

    # Plot results
    # axs[i].scatter(y_pred, y, c=colors)
    axs[i].plot(y_pred, y, "o")

    # Plot parity line
    # axs[i].plot([-4, 10], [-4, 10], "--", color="black")
    axs[i].plot([-4, 10], [-4, 10], "--")

    # Set axis limits
    axs[i].set_xlim(-4, 10)
    axs[i].set_ylim(-4, 10)

    # Add equation
    equation = "$E_v = {:.2f} {:+.2f} E_b {:+.2f} V_r {:+.2f} E_g$".format(cfm.intercept_, cfm.coef_[0], cfm.coef_[1],
                                                                           cfm.coef_[2])
    axs[i].set_xlabel(equation)

    # Add MAE
    mae = mean_absolute_error(y, y_pred)
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
    # axs[i].legend(df_colors['group'].map(color_map), loc='lower right')

    # Add y-axis label
    if i == 0:
        axs[i].set_ylabel("$E_v$ (eV)")

plt.tight_layout()
plt.savefig("kumagai_ternary_vr_from_csv_eb_frac_KF5.png", dpi=300)
plt.show()
