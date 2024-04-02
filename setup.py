from setuptools import setup
__version__ = '0.3'
setup(name='my_routines', author='Joel C. Zinn', install_requires=["astropy", "numpy", "pandas", "scipy", "matplotlib", "emcee", 'corner', 'pudb', 'matplotlib', 'pylab', 'glob', 'mesa_reader'], packages=['./'])
