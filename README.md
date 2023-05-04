# pyoxvac

Pyoxvac is a Python module that defines a Crystal class and some functions that analyze the properties of a crystal structure. 

The Crystal class has three attributes:

    1. `structure`: a pymatgen.Structure object representing the crystal structure.
    2. `unique_oxygens`: a dictionary of unique oxygen sites in the crystal structure, where the keys are the indices of the oxygen sites and the values are pymatgen.Structure objects.
    3. `nonO`: a dictionary of non-oxygen sites in the crystal structure, where the keys are the indices of the non-oxygen sites and the values are pymatgen.Site objects.
    
The class has four methods:

    1. `__init__(self, structure)`: initializes a new Crystal object with the given pymatgen.Structure object.
    2. `visualize(self)`: opens a visualization window of the crystal structure using ASE's view function.
    3. `unique_oxygen(self, structure)`: identifies and stores the unique oxygen sites in the crystal structure as a dictionary of pymatgen.Structure objects.
    4. `non_oxygen(self, structure)': identifies and stores the non-oxygen sites in the crystal structure as a dictionary of pymatgen Site objects.

    
The module also defines four functions:

    1. `file_readin(filepath, visualized)`: reads a crystal structure from a file and returns a Crystal object.
    2. `oxygen_cn(crystal)`: computes the coordination numbers and charges of the unique oxygen atoms in a crystal structure.
    3. `non_oxygen_oxi_state(crystal)`: computes the oxidation states of the non-oxygen atoms in a crystal structure.
    4. `def crystal_data(crystal)`: Generates a dataframe of crystal structure data, in the following format: 
    | Material Name | Index of Unique Oxygen | Coordination Number | Neighbor 1: charge, CN | Neighbor 2: charge, CN | etc. | 
    
    
