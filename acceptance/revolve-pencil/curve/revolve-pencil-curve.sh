#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/egads/one
    two=${HOME}/refine/egads/src
fi

${two}/ref adapt \
    pencil.meshb \
    -g ../gen/pencil.egads \
    -o pencil-curve.meshb \
    -t \
    > accept-annulus-uniform.out

cat accept-annulus-uniform.out
../../check.rb accept-annulus-uniform.out 0.05 3.1

