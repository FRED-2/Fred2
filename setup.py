from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from codecs import open  # To use a consistent encoding
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    readme = f.read()

setup(
    name='FRED2',

    # Version:
    version='2.0.0a1',

    description='A Framework for Epitope Detection',
    long_description=readme,

    # The project's main homepage.
    url='https://github.com/b-schubert/Fred2',

    # Author details
    author='Benjamin Schubert',
    author_email='schubert@informatik.uni-tuebingen.de',

    # Choose your license
    license='BSD',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Biologists, Pharmacologist, Developer',
        'Topic :: Immunoinformatics :: Prediction Tools',

        # The license as you wish (should match "license" above)
        'License :: OSI Approved :: BSD License',

        # The supported Python versions (other than development version were 
        # not yet tested. Especially we should check for Python 3 support
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],

    # What FRED2 relates to:
    keywords='epitope prediction MHC FRED development',

    # Specify  packages via find_packages() and exclude the tests and 
    # documentation:
    packages=find_packages(exclude=['TestSuite', 'doc']),

    # Run-time dependencies. (will be installed by pip when FRED2 is installed)
    install_requires=['pandas' 'coopr' 'biopython'],

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'sample': ['package_data.dat'],
    },

    # 'package_data' is used to als install non package data files
    # see http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
    # example:
    #data_files=[('exapl_target_dir', ['data/file1','data/file2'])],

    # Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'EpitopePrediction=Fred2.Apps.EpitopePrediction:main',
        ],
    },
)
