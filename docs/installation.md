# Installation

## Basic installation

- `pip install dfm_tools`

## Miniforge/conda installation

This is an example with conda, but you can also use any other python environment (plain python, pixi, etc):
- python 3.10 to 3.13 are supported
- download Miniforge3 from [conda-forge.org](https://conda-forge.org/miniforge) and install it with the recommended settings.
- open Miniforge Prompt
- `conda create --name dfm_tools_env python=3.12 git spyder -y` (`git` and `spyder` are optional)
- `conda activate dfm_tools_env`
- install latest dfm_tools release (stable): `pip install dfm_tools -U` (`-U` is for updating)
- to remove environment when necessary: `conda remove -n dfm_tools_env --all`
