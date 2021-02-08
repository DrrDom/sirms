import setuptools
from os import path
import sirms

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name="sirms",
    version=sirms.__version__,
    author="Pavel Polishchuk",
    author_email="pavel_polishchuk@ukr.net",
    description="SiRMS: simplex representation of molecular structure",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DrrDom/sirms",
    packages=['sirms'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    python_requires='>=3.6',
    extras_require={
        'rdkit': ['rdkit>=2017.09'],
    },
    entry_points={'console_scripts':
                  ['sirms = sirms.sirms:main',
                   'calc_atomic_properties_chemaxon = sirms.utilities.calc_atomic_properties_chemaxon:main']}
)
