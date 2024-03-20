from glob import glob
import numpy as np
import pandas as pd
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.core import Structure, Composition
from tqdm import tqdm

from crystal_analysis import Crystal
from utils.data_processing import load_and_merge_csv, filter_defects, add_oxidation_states, prepare_formula_based_features, calculate_crystal_features
from sklearn.linear_model import HuberRegressor

def main():
    data_path = '../data/papers/witman/structures_production/data_01_03_22/'
    csv_pattern = "csvs/*.csv"
    poscar_path = data_path + "poscars"

    # load and merge CSV files
    df = load_and_merge_csv(data_path, csv_pattern, unique_identifier="filename")
    print('step 1 success')

    # filter defects for "V_O"
    df = filter_defects(df, defect_name="V_O")
    print('step 2 success')

    # add oxidation states to structures
    df = add_oxidation_states(df, data_path, oxstate_path_pattern="{filename}_oxstate")
    print('step 3 success')

    # prepare formula-based features
    df = prepare_formula_based_features(df)
    print('step 4 success')

    # calculate crystal features
    df_cf = calculate_crystal_features(df, structure_column="structure")
    print('step 5 success')


    # TODO: refactor code from this point forward to use modeling.py

    # save  cleanedup data to a CSV file
    output_file = "../data/papers/witman/figures/witman_data_cleanup.csv"
    df_cf.to_csv(output_file, index=False)

    # remove NaN values
    df_cf = df_cf.dropna()

    # fit a Huber Regressor model
    cfm = HuberRegressor()
    X = df_cf[["Eb_sum", "Vr_max", "bandgap_eV", "adjusted_dH"]]
    y = df_cf["dH_eV"]
    cfm.fit(X, y)
    y_pred = cfm.predict(X)
    coefs = cfm.coef_

    # print model coefficients and metrics
    print(coefs)
    mae = np.mean(np.abs(y - y_pred))
    print("Mean Absolute Error:", mae)
    print("Sample Size:", len(y))

    # print the equation of the model
    equation = f"$E_v$ = {cfm.intercept_:.2f} + {cfm.coef_[0]:.2f} $\\Sigma E_b$ + {cfm.coef_[1]:.2f} $V_r$ + {cfm.coef_[2]:.2f} $E_g$ + {cfm.coef_[3]:.2f} $E_h_u_l_l$"
    print("Equation:", equation)

if __name__ == "__main__":
    main()