import matplotlib.style
from sklearn.linear_model import HuberRegressor
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from adjustText import adjust_text
from sklearn.model_selection import KFold, train_test_split, cross_val_score
from sklearn.metrics import mean_absolute_error, mean_squared_error

df_cf = pd.read_csv('')
dropped_values = df_cf[df_cf.isnull().any(axis=1)]
dropped_values.to_csv("weighted_drops")

df_cf = df_cf.dropna()
cfm = HuberRegressor()
# X = df_cf[["Vr_max"]]
# X = df_cf[["Eb_sum", "Vr_max"]]
# X = df_cf[["Vr_max", "Eg"]]
X = df_cf[["Eb_sum", "Vr_max", "Eg"]]
# X = df_cf[["Eb_sum", "Vr_max", "Eg", "Ehull"]]
y = df_cf["Ev"]

#cfm
cfm.fit(X, y)
y_pred = cfm.predict(X)
coefs = cfm.coef_

# model tests
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=21)
# # cfm.fit(X_train, y_train)
# # y_pred = cfm.predict(X_test)
# # mae = mean_absolute_error(y_test, y_pred)

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

plt.style.use('seaborn-v0_8-talk')
plt.scatter(y_pred, y)
plt.plot([1, 9], [1, 9], "k--")
plt.xlim(min(y_pred) - 1, max(y_pred) + 1)
plt.ylim(min(y) - 1, max(y) + 1)
plt.gca().set_aspect('equal', adjustable='box')
#equation = f"$E_v$ =  {cfm.coef_[0]:.2f} $V_r$ + {cfm.intercept_:.2f} (eV)"
#equation = f"$E_v$ = {cfm.coef_[0]:.2f} $\\Sigma E_b$ {cfm.coef_[1]:.2f} $V_r$ + {cfm.intercept_:.2f} (eV)"
equation = f"$E_v$ = {cfm.coef_[0]:.2f} $\\Sigma E_b$ {cfm.coef_[1]:.2f} $V_r$ + {cfm.coef_[2]:.2f} $E_g$ + {cfm.intercept_:.2f} (eV)"
# equation = f"$E_v$ = {cfm.coef_[0]:.2f} $\\Sigma E_b$ {cfm.coef_[1]:.2f} $V_r$ + {cfm.coef_[2]:.2f} $E_g$ + {cfm.coef_[3]:.2f} $E_{{hull}}$ + {cfm.intercept_:.2f} (eV)"
mae = np.mean(np.abs(y - y_pred))
plt.text(1.7, 7, f"MAE = {mae:.2f} eV", fontsize=14)
#binary
# oxides = f"$MO_x$"
#binary and ternary
oxides = f"$MO_x$, $ABO_x$"
plt.text(7, 2, oxides, fontsize=14)
#add number of data points as text
plt.text(1.7, 5.5, f"n = {len(y)}", fontsize=14)
plt.text(1.7, 6.5, f"mean score = {-1 * mean:.2f}", fontsize=14)
plt.text(1.7, 6.0, f"std = {std:.2f}", fontsize=14)
texts = []
# for x, y, s in zip(y_pred, y, df_cf["formula"]):
#         texts.append(plt.text(x, y, s, size=10))
# adjust_text(texts, arrowprops=dict(arrowstyle="-", color="k", lw=0.5))
plt.xlabel(str(equation))
plt.ylabel(f"Witman et al. $E_v$")
plt.title("CFM for binary oxides with band gap energies")
#plt.legend(bbox_to_anchor=(1.1, 1.0), prop={'size': 8})
#plt.text(1.1, 0.9, "each color represents a unique defect id")
plt.savefig("CFM_t_Eb_Vr_Eg_nw_KF5.png", dpi=300)
plt.show()