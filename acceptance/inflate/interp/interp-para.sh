#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

np=2

serveCSM -batch poly.csm 

mpiexec -np ${np} ${src}/refmpifull bootstrap poly.egads

mpiexec -np ${np} ${src}/refmpi collar normal \
    poly-vol.meshb \
    10 \
    0.1 \
    2.0 \
    1.68 \
    --fun3d-mapbc poly-vol.mapbc \
    -x surf.meshb

mpiexec -np ${np} ${src}/refmpi \
      adapt poly-vol.meshb \
      -x poly-curve.meshb \
      -f poly-curve-prop.tec

mpiexec -np ${np} ${src}/refmpi collar normal \
    poly-curve.meshb \
    10 \
    0.1 \
    2.0 \
    1.68 \
    --fun3d-mapbc poly-vol.mapbc \
    -x curve.meshb

${src}/ref_acceptance -u u5 surf.meshb surf.solb

mpiexec -np ${np} ${src}/refmpi interp \
	surf.meshb surf.solb curve.meshb curve.solb

exit
