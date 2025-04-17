"""
Support datasets for testing and examples.

Work in progress skeleton functions.
"""

import pooch

DATA_URL = "https://raw.githubusercontent.com/mines-oceanography/"

# file_path = pooch.retrieve(
#     # URL to one of Pooch's test files
#     url="https://github.com/mines-oceanography/ray_tracing/raw/refs/heads/gwen_dev/data/bathy_agulhas.nc",
#     known_hash=None,
# )


def fetch_bathy_agulhas():
    """Download agulhas bathymetry file"""
    file_path = pooch.retrieve(
        # URL to one of Pooch's test files
        url=f"{DATA_URL}swot_wave_current/refs/heads/main/data/Mawar_typhon2.csv",
        known_hash="sha256:647c3809026f6c253af8e079d808dc52ecccd4b333846786895eca93d2268dca",
        progressbar=True,
    )
    return file_path
