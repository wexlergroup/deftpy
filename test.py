import pymatgen.core.structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view
from pymatgen.ext.matproj import MPRester

from mp_api.client import MPRester
with MPRester(api_key="C4VoKwbKXbVC9dk85Uq5G5tig1LWbPKu") as mpr:
    data = mpr.materials.get_data_by_id("mp-4019")

blah = data.structure
blah.add_oxidation_state_by_guess()
print(blah)