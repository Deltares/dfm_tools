# Installation

## Basic installation

- `pip install dfm_tools`

## Miniforge/conda installation

This is an example with conda, but you can also use any other python environment (plain Python venv, pixi, uv etc):
- python 3.11 to 3.14 are supported
- download Miniforge from https://conda-forge.org/download and install it with the recommended settings.
- open Miniforge Prompt
- `conda create --name dfm_tools_env python=3.13 -y`
- `conda activate dfm_tools_env`
- install latest dfm_tools release (stable): `pip install dfm_tools -U` (`-U` is for updating)
- to remove environment when necessary: `conda remove -n dfm_tools_env --all`
