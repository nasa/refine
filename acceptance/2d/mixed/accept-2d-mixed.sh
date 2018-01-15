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

ref_acceptance square_0.msh square_0.metric 0.01
ref_driver -i square_0.msh -m square_0.metric -o mixed01 > accept-2d-mixed-01.status

ref_acceptance mixed01.b8.ugrid mixed01.metric 0.01
ref_driver -i mixed01.b8.ugrid -m mixed01.metric -o mixed02 > accept-2d-mixed-02.status

cat accept-2d-mixed-02.status
../../../check.rb accept-2d-mixed-02.status 0.60 1.5

