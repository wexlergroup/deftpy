import pymatgen.core.structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher
from ase.visualize import view
from pymatgen.ext.matproj import MPRester

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
def blah(*args):
    structure = pymatgen.core.Structure.from_file(args[0])
    structure_matcher = pymatgen.analysis.structure_matcher.StructureMatcher()

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
    unique_structures[0] = vac_structures[0]  
    
    for i in vac_structures[1:]:
        n_dupl = 0
        for j in unique_structures[1:]:
            if structure_matcher.fit(vac_structures[i], vac_structures[j]) == True:
                n_dupl += 1
        if n_dupl == 0:
            unique_structures[i] = vac_structures[i]
        
    return(unique_structures)
    
data = blah("mk1/crystal_files/OQMD_CaTiO3_POSCAR.txt")
print(data)







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