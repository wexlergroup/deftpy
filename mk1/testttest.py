import pymatgen.core.structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher
from ase.visualize import view
from pymatgen.ext.matproj import MPRester
#import CrystalNN
from pymatgen.analysis.local_env import CrystalNN
import pandas as pd
from pymatgen.core.composition import Composition


structure = pymatgen.core.Structure.from_file('mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt')
structure_matcher = pymatgen.analysis.structure_matcher.StructureMatcher()
structure.add_oxidation_state_by_guess()

iO = []
vac_structures = {}
for i, site in enumerate(structure.sites):
    if site.specie.symbol == 'O':
        iO.append(i)
for i in iO:
    new_struct = structure.copy()
    new_struct.remove_sites([i])
    vac_structures[i] = new_struct

unique_structures = {}

for i, vac_structure in vac_structures.items():
    if i == 0:
        unique_structures[i] = vac_structure
    n_dupl = 0
    for j, unique_structure in unique_structures.items():
        if structure_matcher.fit(vac_structure, unique_structure) == True:
            n_dupl += 1
    if n_dupl == 0:
        unique_structures[i] = vac_structure


struc_NN = {}

for i in range(len(unique_structures.keys())):
    pos_arg = list(unique_structures.keys())[i]

    struc_NN[i] = CrystalNN(weighted_cn=True, cation_anion=True, porous_adjustment=False, \
    distance_cutoffs=None, x_diff_weight=0).get_nn_info(structure, pos_arg)



print(struc_NN.values())

df_struc = pd.DataFrame(struc_NN)
print(df_struc.head())

pd.DataFrame.to_csv(df_struc, 'mk1/crystal_files/df_struc.csv')

for i in df_struc:
    