#!/usr/bin/python
"""
WARNING: GEOPACK must be compiled with 8-byte reals!  See README for
more information.
"""
import sys
from numpy.distutils.core import setup, Extension

geopack = Extension(name = 'geopack', sources = ['geopack.pyf','geopack.F'])

if __name__ == "__main__":
    setup(name = 'geopack',
          version = '1.0',
          description       = 'Python bindings for GEOPACK 2005',
          author            = 'Peter Schmitt',
          author_email      = 'schmitt@ucar.edu',
          ext_modules = [geopack]
          )
