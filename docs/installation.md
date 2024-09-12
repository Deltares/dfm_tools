# Installation

## Basic installation

- `pip install dfm_tools`

## Recommended installation

- download Miniforge3 from [the miniforge github](https://github.com/conda-forge/miniforge?tab=readme-ov-file#download) and install it with the recommended settings.
- open Miniforge Prompt
- `conda create --name dfm_tools_env python=3.12 git spyder -y` (`git` and `spyder` are optional)
- `conda activate dfm_tools_env`
- install latest dfm_tools release (stable): `pip install dfm_tools -U` (`-U` is for updating)
- alternatively install most recent dfm_tools version from github (might be unstable): `python -m pip install git+https://github.com/Deltares/dfm_tools`
- to remove environment when necessary: `conda remove -n dfm_tools_env --all`

## Potential git/conda issues with installation

- long paths error? Check [this Github issue](https://github.com/Deltares/HYDROLIB-core/issues/327#issuecomment-1266534032)
- OpenSSL error? Fix your conda base env by doing [this](https://github.com/conda/conda/issues/11795#issuecomment-1335666474) or maybe [this](https://github.com/conda/conda/issues/11795#issuecomment-1382661765). Let us know if you encounter this issue.
