from setuptools import setup
__version__ = '0.2'
setup(name='my_routines', author='Joel C. Zinn', install_requires=["pyfits", "numpy", "pandas", "scipy", "matplotlib", "emcee", 'corner', 'pudb'], packages=['./'])
