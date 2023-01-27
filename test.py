from mp_api.client import MPRester
with MPRester(api_key="C4VoKwbKXbVC9dk85Uq5G5tig1LWbPKu") as mpr:
    data = mpr.materials.get_data_by_id("mp-4019")

print(data)