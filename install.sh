#!/bin/sh

# Installs dependencies required by pyLTR
# Warning: Only tested with standard Linux environments!

cd pkgs/cxform-0.71_python
python setup-cxform.py install --install-lib=../../ghostpy/transform
cd ../..

### 
### # Note: GEOPACK must be compiled with 8-byte reals.  The flag to do
### # this is compiler dependent!  See pkgs/geopack/README for details.
### cd pkgs/geopack-2005
### python setup-geopack.py config_fc --fcompiler=gnu95 --f77flags=-fdefault-real-8 install --install-lib=../../pyLTR/transform
### cd ../..

cd pkgs/geopack-2008
python setup-geopack.py config_fc --fcompiler=gnu95 --f77flags=-fdefault-real-8 install --install-lib=../../ghostpy/transform
cd ../..
