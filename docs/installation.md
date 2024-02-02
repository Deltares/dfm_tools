# Installation

## Basic installation

- ``pip install dfm_tools``

## Recommended installation

- python 3.11 is recommended, python>=3.9 is required ([more info](https://github.com/Deltares/dfm_tools/issues/267)), python 3.13 not yet supported ([more info](https://github.com/Deltares/dfm_tools/issues/559))
- download and install Anaconda from [anaconda.com](https://www.anaconda.com)
- open Anaconda prompt
- ``conda create --name dfm_tools_env python=3.11 git spyder -c conda-forge -y`` (``git`` and ``spyder`` are optional)
- ``conda activate dfm_tools_env``
- install latest dfm_tools release (stable): ``pip install dfm_tools -U`` (the ``-U`` is for updating)
- alternatively install most recent dfm_tools version from github (might be unstable): ``python -m pip install git+https://github.com/Deltares/dfm_tools``
- to remove environment when necessary: ``conda remove -n dfm_tools_env --all``

## Potential git/conda issues with installation

- long paths error? Check [this Github issue](https://github.com/Deltares/HYDROLIB-core/issues/327#issuecomment-1266534032)
- OpenSSL error? Fix your conda base env by doing [this](https://github.com/conda/conda/issues/11795#issuecomment-1335666474) or maybe [this](https://github.com/conda/conda/issues/11795#issuecomment-1382661765). Let us know if you encounter this issue.
