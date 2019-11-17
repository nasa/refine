#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/strict/one
    two=${HOME}/refine/strict/src
fi

metric="-twod linear-0001"

${two}/ref_acceptance ${metric} square_0.msh square_0.solb
${two}/ref_driver -i square_0.msh -m square_0.solb -x mixed01.meshb -t > accept-2d-mixed-01.status

${two}/ref_acceptance ${metric} mixed01.meshb mixed01.solb
${two}/ref_driver -i mixed01.meshb -m mixed01.solb -x mixed02.meshb -t > accept-2d-mixed-02.status

cat accept-2d-mixed-02.status
../../check.rb accept-2d-mixed-02.status 0.1 12.0

