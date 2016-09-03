#!/usr/bin/env python

# this is a generic installation file for a python code, taken from Daniel Foreman-Mackey.

import re
import os
import sys
from setuptools import setup

module = "my_routines"

version = re.findall(r"__version__ = \"(.*?)\"", open(module+".py").read())[0]

setup(
    name=module,
    version=version,
    author="Joel C. Zinn",
    author_email="zinn.44@osu.edu",
#    url="https://github.com/dfm/triangle.py",
    py_modules=[module],
#    description="",
#    long_description=open("README").read(),
#    package_data={"": ["LICENSE"]},
#    include_package_data=True,
#    classifiers=[
#        "Development Status :: 5 - Production/Stable",
#        "License :: OSI Approved :: BSD License",
#        "Intended Audience :: Developers",
#        "Intended Audience :: Science/Research",
#        "Operating System :: OS Independent",
#        "Programming Language :: Python",
#    ],
)

