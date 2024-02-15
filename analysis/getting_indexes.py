import pandas as pd
import os

#path for scratch slightly different than for analysis scripts
data_path = "/Users/isakov/PycharmProjects/deft_ethans_branch/deftpy/data/papers/kumagai/"

# Get data
df_0 = pd.read_csv(data_path + "vacancy_formation_energy_ml/charge0.csv")  # neutral vacancies
df_1 = pd.read_csv(data_path + "vacancy_formation_energy_ml/charge1.csv")  # +1 charged vacancies
df_2 = pd.read_csv(data_path + "vacancy_formation_energy_ml/charge2.csv")  # +2 charged vacancies

# Add charge column
df_0["charge"] = 0
df_1["charge"] = 1
df_2["charge"] = 2

# Combine dataframes
df = pd.concat([df_0, df_1, df_2], ignore_index=True).reset_index(drop=True)

# Remove the column named "Unnamed: 0"
df = df.drop("Unnamed: 0", axis=1)

# Retrieve index sites for the vacancy in the supercell.vesta file
vacancy_indexes = []
for defect in (df["vacancy_formation_energy"].unique()):
    df_defect = df[df["vacancy_formation_energy"] == defect]
    full_name = df_defect["full_name"].iloc[0]
    formula = full_name.split("_")[0]
    vacancy = full_name.split("Va_")[1]
    info_path = os.path.join(data_path, "site_info", formula, "cell_info.txt")
    with open(info_path, "r") as txt:
        index = "Irreducible element: " + vacancy
        lines = txt.readlines()
        for i, line in enumerate(lines):
            if index in line:
                target_line_index = min(i + 5, len(lines) - 1)
                target_line = lines[target_line_index]
                vacancy_index = int(target_line.split("..")[-1])
                vacancy_index = vacancy_index
                vacancy_indexes.append(vacancy_index)
                break

# Add vacancy index data to the complete df
df["vacancy_index"] = vacancy_indexes
print(df.columns)

df.to_csv("complete_df_indexed.csv")