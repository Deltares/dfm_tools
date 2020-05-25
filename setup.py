#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="Jelmer Veenstra",
    author_email='jelmer.veenstra@deltares.nl',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: Beta',
        'Intended Audience :: Modellers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="dfm_tools are post-processing tools for Delft3D FM",
    entry_points={
        'console_scripts': [
            'dfm_tools=dfm_tools.cli:main',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='dfm_tools',
    name='dfm_tools',
    packages=find_packages(include=['dfm_tools', 'dfm_tools.io']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/openearth/dfm_tools',
    version='0.7.21',
    zip_safe=False,
)
