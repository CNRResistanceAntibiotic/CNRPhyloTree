#!/usr/bin/python3

# setup.py file install libraries wanted for cnrphylotree
# launch it with pip install -e .

# Install setuptools if not already present.
from setuptools import setup, find_packages
import glob

setup(
    name='cnrphylotree',
    version='1.0.1',
    description='cnrphylotree: pipeline CNR Resistance for phylo tree',
    packages=find_packages(),
    package_data={
        'cnr_phylo_tree_src': ['save_newick_to_png.R'],
    },
    author='Richard Bonnet',
    author_email='rbonnet@chu-clermontferrand.fr',
    url='https://github.com/CNRResistanceAntibiotic/CNRPhyloTree',
    scripts=glob.glob('scripts/*'),
    install_requires=[''],
    license='GPLv3',
    classifiers=[
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],

)
