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

serveCSM -batch poly.csm 

${two}/ref bootstrap poly.egads

rm -rf inflated.b8.ugrid
${two}/ref_inflatable  \
    poly-vol.meshb \
    10 \
    0.1 \
    2.0 \
    1.68 \
    2 3 4 5 6 7 8 9

${two}/ref translate inflated.b8.ugrid surf.meshb

${two}/ref \
      adapt poly-vol.meshb \
      -g poly.egads \
      -x poly-curve.meshb \
      -f poly-curve-prop.tec

rm -rf inflated.b8.ugrid
${two}/ref_inflatable  \
    poly-curve.meshb \
    10 \
    0.1 \
    2.0 \
    1.68 \
    2 3 4 5 6 7 8 9

${two}/ref translate inflated.b8.ugrid curve.meshb

${two}/ref_acceptance -u 5 surf.meshb surf.solb

${two}/ref interp surf.meshb surf.solb curve.meshb curve.solb

exit
