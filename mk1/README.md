The goal for mk1 is to read in a crystal structure file (CSF) from various databases (MP, OQMD, AFLOW, ICSD, etc). into pymatgen. From pymatgen, we will use the pymatgen.io.ase module to pipeline the CSF to Atomic Simulation Environment. Using the view module from ASE, we will print the crystal structure and test it against a printed VESTA crystal structure of the same variety.

Working backwards, the first test case we will aim to fulfill is to match the ASE CSF with the VESTA CSF.

Further tests we will employ will be within the structure itself, such as the number of neighbors to an oxygen or the specific neighbors to an oxygen.

Packages used: pymatgen, ASE, VESTA
