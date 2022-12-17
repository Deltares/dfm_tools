#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    author="Jelmer Veenstra",
    author_email='jelmer.veenstra@deltares.nl',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: Beta',
        'Intended Audience :: Modellers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="dfm_tools are pre- and post-processing tools for Delft3D FM",
    entry_points={
        'console_scripts': [
            'dfm_tools=dfm_tools.cli:main',
        ],
    },
    install_requires=['scipy', 'numpy', 'matplotlib', 'pandas', 'hydrolib-core>=0.3.1', 'shapely', 'cartopy', 'pyepsg', 'geopandas', 'contextily', 'xarray', 'dask', 'netcdf4>=1.5.3', 'bottleneck', 'xugrid', 'cdsapi', 'pydap'],
    license="GNU General Public License v3",
    long_description=readme,
    include_package_data=True,
    keywords='dfm_tools',
    name='dfm_tools',
    packages=find_packages(include=['dfm_tools', 'dfm_tools.io']),
    test_suite='tests',
    tests_require=['pytest>=5.0.1','bump2version>=0.5.11','pytest-cov','pdoc3','pandoc'],
    url='https://github.com/Deltares/dfm_tools',
    version='0.8.36',
    zip_safe=False,
)
