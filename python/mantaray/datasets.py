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


Maverick = pooch.create(
    path=pooch.os_cache("mantaray"),
    base_url="https://github.com/mines-oceanography/ray_tracing/releases/download/d0.0.7/",
    registry={
        "bathy_agulhas.nc": "sha256:706500b740f5212176b7598e8f7dddabcff3d9d2c21eae744a492236b8ba2df8",
        "bathy_nazare_edited_bis.nc": "ha256:d48a47263f90a8fd40c2ea4e51ed289b1c5bf7b57ae410fdbfff0e97cb5d9402",
        "current_agulhas.nc": "sha256:13bad7ec6be48a1ebf0c11e1f1a5f3c2ce0bbb5110598276d9859117c80aeec6",
        "current_nazare.nc": "sha256:4471d887fbee0a851f060bbd95342a57194b9dae4f17434ea4b1ce3aef469c84",
        "globcurrent.nc": "sha256:8631c7c8f74f7d460557105f79daf86a8130ce4617747aac159be9743c48c194",
        "nazare_bathy.nc": "sha256:ca38d35f2a686de3daf1865997c5817bf5bc6189fdf5b752ac997fbe35b1b7b0",
        "s2_image.nc": "sha256:95d202ac05361c14356273ab323ea5b3e7bbb6a11e09e28db7807c460caa7bba",
        "spec_s2.nc": "sha256:c5463822f7da499b99d5e870cbb30de263835e16cba2380e8148c26ab2aab519",
        "wavewatch_bathy.nc": "sha256:3a2292171626a94aaf7ed05492bf97d4884bef6967425c2629da22521f4ab8ad",
        "wavewatch_currents.nc": "sha256:dfcd89d24f52d71c24853ce1d77f5e67cccc81bdda384bbd139f3bb5a7e9e55b",
    },
)

def fetch_bathy_agulhas():
    """Download agulhas bathymetry file"""
    file_path = pooch.retrieve(
        # URL to one of Pooch's test files
        url=f"{DATA_URL}swot_wave_current/refs/heads/main/data/Mawar_typhon2.csv",
        known_hash="sha256:647c3809026f6c253af8e079d808dc52ecccd4b333846786895eca93d2268dca",
        progressbar=True,
    )
    return file_path
