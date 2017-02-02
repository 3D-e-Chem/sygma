import os
from setuptools import setup, find_packages

# from setuptools.extension import Extension

here = os.path.abspath(os.path.dirname(__file__))
try:
    README = open(os.path.join(here, 'README.rst')).read()
except IOError:
    README = ''

# Test if rdkit is present with INCHI support
try:
    from rdkit.Chem.inchi import INCHI_AVAILABLE

    if not INCHI_AVAILABLE:
        raise Exception('RDKit with INCHI support is required')
except ImportError:
    raise Exception('RDKit with INCHI support is required')

ext_modules = []

setup(
    name='SyGMa',
    version="1.0.1",
    license='GPL',
    author='Lars Ridder',
    author_email='lars.ridder@esciencecenter.nl>',
    url='https://github.com/ridderl/sygma',
    description='Systematic Generation of potential MetAbolites',
    long_description=README,
    classifiers=["Intended Audience :: Science/Research",
                 "Environment :: Console",
                 "Natural Language :: English",
                 "Operating System :: OS Independent",
                 "Topic :: Scientific/Engineering :: Chemistry",
                 "Programming Language :: Python :: 2.7",
                 ],
    packages=find_packages(),
    package_data={'sygma': ['rules/*.txt']},
    entry_points={'console_scripts': ['sygma = sygma.script:main']}
)
