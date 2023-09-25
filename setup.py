from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("generate_kmers_dna.pyx")
)
