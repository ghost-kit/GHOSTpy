#! /bin/bash

fileName=$1
vector=$2
radius=$3

python /glade/u/home/jmurphy/src/ghostpy/samples/parallel_lstar/parallelLstar.py -f $fileName  -v $vector -r $radius

