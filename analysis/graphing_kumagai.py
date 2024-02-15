import pandas as pd

import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.metrics import mean_absolute_error

df_plot = pd.read_csv("kumagai_Eb_Vr.csv")

fig, axs = plt.subplots(ncols=3, figsize=(12, 4))

for i, charge in enumerate([0, 1, 2]):
    cfm = linear_model.HuberRegressor()
    #X = df_plot.loc[df_plot.charge == charge, ["vr_max", "band_gap"]]
    #y = df_plot.loc[df_plot.charge == charge, "vacancy_formation_energy"]
    X = df_plot.loc[df_plot.charge == charge, ["Eb_sum", "vr_max", "band_gap"]]
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
    equation = "$E_v = {:.2f} {:+.2f} E_b {:+.2f} V_r {:+.2f} E_g$".format(cfm.intercept_, cfm.coef_[0], cfm.coef_[1],
                                                                           cfm.coef_[2])
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
plt.savefig("kumagai_ternary_vr_from_csv_eb.png", dpi=300)
plt.show()