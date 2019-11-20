#! /bin/bash

fileName=$1
vector=$2
radius=$3

mpirun -np 24 python ./parallelLstar.py -f $fileName  -v $vector -r $radius

