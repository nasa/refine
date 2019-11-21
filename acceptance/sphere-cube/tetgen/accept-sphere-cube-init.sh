#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

cp ../geom/sphere-cube.egads .

${src}/ref boostrap sphere-cube.egads \
      | tee sphere-cube-init.out

${src}/ref adapt sphere-cube-vol.meshb -g sphere-cube.egads \
      -x sphere-cube-crv.meshb -t | tee sphere-cube-crv.out

