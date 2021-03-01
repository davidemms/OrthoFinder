# HOW TO SETUP (FOR PIP)
# Go to the orthofinder directory and remove dist/ and build folders by
# sudo rm -r dist/
# sudo rm -r build/
# Create the bdist_wheel file:
# sudo python setup.py sdist bdist_wheel
# Upload to pypi:
# twine upload dist/*
# You will be asked for username (davidemms) and password
# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
import sys
from os import path
from scripts_of.util import version


setup(
    name='orthofinder',
    version=version,

    description='Phylogenetic orthology inference for comparative genomics',
    long_description='Phylogenetic orthology inference for comparative genomics',

    # The project's main homepage.
    url='https://github.com/davidemms/OrthoFinder',

    # Author details
    author='David Emms',
    author_email='david_emms@hotmail.com',

    # Choose your license
    license='OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',

    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',

        # Specify the Python versions you support here.
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    # What does your project relate to?
    keywords='',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['docs', 'tests']),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['numpy', 'scipy'],

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
     package_data={
          '': ['config.json'],
        #   'ExampleData': ['ExampleData/Mycoplasma_hyopneumoniae.faa', 'ExampleData/Mycoplasma_agalactiae.faa', 'ExampleData/Mycoplasma_gallisepticum.faa', 'ExampleData/Mycoplasma_genitalium.faa']
     },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    data_files=[('ExampleData', ['ExampleData/Mycoplasma_hyopneumoniae.faa', 'ExampleData/Mycoplasma_agalactiae.faa', 'ExampleData/Mycoplasma_gallisepticum.faa', 'ExampleData/Mycoplasma_genitalium.faa'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'orthofinder=scripts_of.__main__:main',
            'primary_transcript=tools.primary_transcript:main',
            'make_ultrametric=tools.make_ultrametric:main',
            'convert_orthofinder_tree_ids=tools.convert_orthofinder_tree_ids:main'
        ],
    },
)
