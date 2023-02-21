import pymatgen.core.structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher
from ase.visualize import view
from pymatgen.ext.matproj import MPRester
#import CrystalNN
from pymatgen.analysis.local_env import CrystalNN
import pandas as pd
from pymatgen.core.composition import Composition
# for using data from the Materials Project, enter a tuple with the 
# mp-id as the first element and your api key as the second element

print('hi')
def file_readin(*args):
    if args[0].startswith('mp-'):
        with MPRester(api_key=(args[1])) as mpr:
            data = mpr.materials.get_data_by_id(args[0]) 
    else:
        structure = pymatgen.core.Structure.from_file(args[0])
        print( structure.sites())
        visualized = AseAtomsAdaptor.get_atoms(structure)
        #works with cif, vasp, and poscar files

    return(view(visualized))


print('hi2')
def oxygenCN(*args):
    structure = pymatgen.core.Structure.from_file(args[0])
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
        # new_struct.add_oxidation_state_by_guess()
        vac_structures[i] = new_struct

    unique_structures = {}
    
    # unique_structures[0] = vac_structures[0]  
    print(vac_structures[4])
    for i, vac_structure in vac_structures.items():
        if i == 0:
            unique_structures[i] = vac_structure
        n_dupl = 0
        for j, unique_structure in unique_structures.items():
             if structure_matcher.fit(vac_structure, unique_structure) == True:
                n_dupl += 1
        if n_dupl == 0:
            unique_structures[i] = vac_structure

# use this


    # for i in vac_structures[1:]:
    #     n_dupl = 0
    #     for j in unique_structures[1:]:
    #         if structure_matcher.fit(vac_structures[i], vac_structures[j]) == True:
    #             n_dupl += 1
    #     if n_dupl == 0:
    #         unique_structures[i] = vac_structures[i]
        
        struc_NN = {}
        for i in range(len(unique_structures.keys())):
            pos_arg = list(unique_structures.keys())[i]
            
            struc_NN[i] = CrystalNN(weighted_cn = True, cation_anion = True, porous_adjustment = False, distance_cutoffs = None, x_diff_weight = 0).get_nn_info(structure, pos_arg)
        

    crystal_df = pd.DataFrame.from_dict(struc_NN, orient='index')
    #   get coordination number subset by the species its bonded to
    
    crystal_df.to_csv('mk1/crystal_files/test.csv', index=False)

        # speed  - 359
    return(crystal_df)
data = oxygenCN('mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt')
print(data)

    
# data = oxygenCN("mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt")
# print(data)

# get the unique periodic sites



# for each index of oxygen thats unique, he wants to see something about the number of neighbor types 
# (in CaTi)3 thatd be 2 (ca/ti) and a list of neighbor 1, its name

#material name, ox4,   # of neighbor types (2) , neighbor 1 , charge of n1, CN 1, neighbor 2, charge of n2, CN 2, etc

# get this and put it into a dataframe.


# every time i add something to the dataframe, i dont want to have to do it again, maybe storing 
# it in the df and checking if its there before doing it again

CeMnO3

def non_o_oxidation_state(*args):
    structure = pymatgen.core.Structure.from_file(args[0])
    notO = []
    for i, site in enumerate(structure.sites):
        if site.specie.symbol != 'O':
            notO.append(i) 
    
    oxi_states = {}
    for j in notO:
        oxi_states[j] = structure[j].specie.add_charges_from_oxi_state_guesses

    return(oxi_states)
# non_o_oxidation_state("mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt")

# add_charges_from_oxi_state_guesses
    


# structure = pymatgen.core.Structure.from_file("mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt")
# wu = []
# for site in (structure.sites):
#     if site.specie.symbol == 'O':
#         wu.append(site)
# oxi_states = {}
# for j in wu:
#     oxi_states[j] = structure[j].specie.add_charges_from_oxi_state_guesses
# print(oxi_states)



        # unique_structures[0] = struc without 0 in dic

    # for i in not 0:
    #     n_dupl = 0
    #     for j in unique_structures:
    #         if structure_matcher.fit(new, old) == True:
    #             n_dupl += 1
    #     if n_dupl == 0:
    #         unique_structures[i] = vac_struct.values()[i]

            # or nested for


        # loop over structures in unique_structures:
        #     if structure_matcher.fit(new, old) != True:
        #         unique_structures[i] = new_struct[i]
        #     else:
        #         unique_structures[i] = vac_structures[i]
    



    # vac_list = []
    # for i in vac_structures:
    #     vac_list.append(structure_matcher.group_structures(vac_structures[i], structure))
    # return(vac_list)