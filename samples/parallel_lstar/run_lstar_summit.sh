#! /bin/bash

fileName=$1
vector=$2
radius=$3

python /projects/jomu9721/src/ghostpy/samples/parallel_lstar/parallelLstar.py -f $fileName  -v $vector -r $radius

