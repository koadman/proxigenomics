#!/bin/bash

#
# Batch script for doing all scores from one invocation in scons.
#

bin/f1score.py $1 $2 ${2}.f1
bin/vmeasure.py $1 $2 ${2}.vm
bin/bcubed.py -o ${2}.bc $1 $2
