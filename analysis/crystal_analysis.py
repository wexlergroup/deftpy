import re
import sys
from enum import Enum
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Callable

import numpy as np
import pandas as pd
from ase.visualize import view
from pymatgen.analysis.defects.generators import VacancyGenerator
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.core import Species, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.bond_valence import BVAnalyzer as BVA
from pymatgen.io.vasp import Poscar

EB_DICT = {"filepath": "../data/features/Eb.csv", "column_name": "Eb", "comparison": "os"}
VR_DICT = {"filepath": "../data/features/Vr.csv", "column_name": "Vr", "comparison": "n"}


class Crystal:
    """
    A class for analyzing crystal structures.

    Attributes:
        structure: The pymatgen Structure object.
        nn_finder: The CrystalNN object.
        eb: The bond dissociation enthalpies.
        vr: The reduction potentials.
        cn_dicts: A list of coordination number dictionaries.
        bond_dissociation_enthalpies: A list of bond dissociation enthalpies.
        reduction_potentials: A list of reduction potentials.


    Methods:
        visualize: Visualizes the crystal structure using ASE's view function.


    Examples:
        # TODO: Add examples
    """

    @staticmethod
    def _split_before_first_number(s: str) -> List[str]:
        """
        Splits a string before the first number.

        Args:
            s: The string to split.

        Returns:
            A list of strings.

        Examples:
            # TODO: Add examples
        """
        return re.split(r"(?=\d)", s, maxsplit=1)

    @staticmethod
    def _parse_species_string(species_string: str) -> Tuple[Optional[Species], str, int]:
        """
        Parses a species string.

        Args:
            species_string: The species string to parse.

        Returns:
            A tuple of the species, the symbol, and the oxidation state.

        Examples:
            # TODO: Add examples
        """
        if species_string == "":
            return None, "", 0
        
        # Check if the string is of valid species format
        if re.match(r"[A-Za-z]+\d+\+", species_string):
            species = Species.from_str(species_string)
            return species, species.symbol, species.oxi_state
        else:
            # Handle strings without numbers or not in expected format
            split_str = Crystal._split_before_first_number(species_string)
            symbol = split_str[0] if split_str else species_string
            oxi_state = round(float(split_str[1][:-1])) if len(split_str) > 1 else 0
            return None, symbol, oxi_state

    def __init__(
            self,
            filepath: Optional[str] = None,
            poscar_string: Optional[str] = None,
            pymatgen_structure: Optional[Structure] = None,
            nn_finder: Optional[CrystalNN] = None,
            # bv_finder: Optional[BVA] = None,
            use_weights: Optional[bool] = False,
            species_symbol: Optional[str] = "O",
            n: Optional[int] = None                  # attempt to allow for pre-emptive indexing
    ):
        
        """
        Initializes the Crystal object.

        Args:
            filepath: The filepath to the POSCAR file.
            poscar_string: The POSCAR string.
            pymatgen_structure: The pymatgen Structure object.
            nn_finder: The CrystalNN object.

        Raises:
            ValueError: If neither filepath, poscar_string, nor pymatgen_structure is specified.

        Examples:
            # TODO: Add examples
        """
        if filepath:
            self.structure = Structure.from_file(filepath)
            if filepath.endswith(".txt"):  # Assuming your POSCAR files have a .txt extension
                self.structure.add_oxidation_state_by_guess()
        elif poscar_string:
            self.structure = Structure.from_str(poscar_string, fmt="poscar")
            self.structure.add_oxidation_state_by_guess()
        elif pymatgen_structure:
            self.structure = pymatgen_structure
        else:
            raise ValueError("Specify either filepath, poscar_string, or pymatgen_structure.")

        self.nn_finder = nn_finder or CrystalNN()
        self.use_weights = use_weights
        self.n = n
        # self.bv_finder = bv_finder
        self.species_symbol = species_symbol
        package_dir = Path(__file__).parent
        self.eb = pd.read_csv(package_dir / EB_DICT["filepath"])
        self.vr = pd.read_csv(package_dir / VR_DICT["filepath"])
        self._cn_dicts_initialized = False
        self.cn_dicts = []
        # self.bv_dicts = []

        self.bond_dissociation_enthalpies = self._get_values(self.eb, EB_DICT["column_name"], EB_DICT["comparison"])


        self.reduction_potentials = self._get_values(self.vr, VR_DICT["column_name"], VR_DICT["comparison"])


    def _initialize_structure_analysis(self) -> List[Dict[str, int]]:
        """
        Initializes the structure analysis.

        Returns:
            A list of coordination number dictionaries.

        Examples:
            # TODO: Add examples
        """
        if self._cn_dicts_initialized:
            return self.cn_dicts

        # Check for oxidation states and add them if they are not present in the structure object already
        if sum([x.oxi_state != 0 for x in self.structure.species]) == 0:
            self.structure.add_oxidation_state_by_guess()
        vacancy_generator = VacancyGenerator()
        # if available, take in a vacancy instead of generating (bottleneck)
        if self.n is not None:
            self.cn_dicts = [self.nn_finder.get_cn_dict(self.structure, n=self.n, use_weights=self.use_weights)]
            # self.cn_dicts = self.nn_finder.BVA().get_valences(self.structure) # DOESNT WORK
            self._cn_dicts_initialized = True
            return self.cn_dicts
            # self.bv_dicts = [self.bv_finder.get_valences(structure=self.structure)]
        else:
        # bulk visual data json
            vacancies = vacancy_generator.get_defects(self.structure)
            indices = [v.defect_site_index for v in vacancies if v.site.specie.symbol == self.species_symbol]
            self.cn_dicts = [self.nn_finder.get_cn_dict(self.structure, i, use_weights=self.use_weights) for i in indices]
            self._cn_dicts_initialized = True
            return self.cn_dicts

    def _get_values(self, dataframe: pd.DataFrame, column_name: str, comparison: str) -> List[Dict[str, float]]:
        """
        Gets the values from the dataframe.

        Args:
            dataframe: The dataframe.
            column_name: The column name.
            comparison: The comparison column name.

        Returns:
            A list of dictionaries.

        Examples:
            # TODO: Add examples
        """
        try:
            self._initialize_structure_analysis()
            values = []
            for cn_dict in self.cn_dicts:
                value = {}
                for species_string, cn in cn_dict.items():
                    species = Species.from_str(species_string)
                    symbol = species.symbol
                    oxidation_state = species.oxi_state
                    condition = (dataframe.elem == symbol) & (dataframe[comparison] == oxidation_state)
                    if not dataframe.loc[condition, column_name].empty:
                        value[species_string] = dataframe.loc[condition, column_name].iloc[0]
                    else:
                        value[species_string] = np.nan
                values.append(value)
            return values
        except KeyError as e:
            raise ValueError(f"Missing required column in dataframe: {e}")
        except IndexError:
            raise ValueError("Index error occurred while accessing dataframe")
        except Exception as e:
            raise ValueError(f"An unexpected error occurred: {e}")
        return values


    def visualize(self):
        """
        Visualizes the crystal structure using ASE's view function.

        Examples:
            # TODO: Add examples
        """
        atoms = AseAtomsAdaptor.get_atoms(self.structure)
        view(atoms)

    def __repr__(self) -> str:
        """
        Returns the string representation of the Crystal object.

        Returns:
            The string representation of the Crystal object.

        Examples:
            # TODO: Add examples
        """
        return f"Crystal({self.structure})"