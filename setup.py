import sys
from setuptools import setup, find_packages
from Cython.Build import cythonize
from PanPA.version import __version__

'''
I can also check for cython's version using
from Cython.Compiler.Version import version

now version is a string, e.g. 0.29.21
'''

'''
It is better to ship this without the need to have Cython
I should do what is described here
https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html
at section Distributing Cython Modules

The idea is that I use extensions and have a check whether I am requiring cython or not
'''

CURRENT_PYTHON = sys.version_info[:2]
REQUIRED_PYTHON_LOWER = (3, 6)
REQUIRED_PYTHON_UPPER = (3, 8)

if CURRENT_PYTHON < REQUIRED_PYTHON_LOWER or CURRENT_PYTHON > REQUIRED_PYTHON_UPPER:
    sys.stderr.write("PanPA requires Python betwee 3.6 and 3.8 "
                     "you current verions is {}".format(CURRENT_PYTHON))
    sys.exit(1)

with open("README.md", "r") as readme:
    long_description = readme.read()

reqs = []
with open("requirements.txt", "r") as infile:
    for l in infile:
        if not l.startswith("#"):
            reqs.append(l.strip())

setup(
    name="PanPA",
    version=__version__,
    license="MIT",
    author="Fawaz Dabbaghie",
    url='https://fawaz-dabbaghieh.github.io/',
    description="Building and aligning to protein graphs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # keywords="proteins alignment graphs pangenome bioinformatics software",
    classifiers=[
        "Development Status :: 3 - Alpha",
        # "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Programming Language :: Python :: 3.6",
        "Operating System :: POSIX :: Linux",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bioinformatics",
    ],
    setup_requires=[],
    tests_require=['pytest'],
    include_package_data=True,
    python_requires=">=3.6",
    packages=find_packages(),
    install_requires=[],
    ext_modules=cythonize("PanPA/*pyx", compiler_directives={"boundscheck": False, "cdivision": True,
                                               "nonecheck": False, "initializedcheck": False,
                                               "language_level": "3"}),
    entry_points={
        "console_scripts": ["PanPA = PanPA.main:main"],
    },
)
