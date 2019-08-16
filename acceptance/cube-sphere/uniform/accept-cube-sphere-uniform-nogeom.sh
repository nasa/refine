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

${two}/ref_translate cube-sphere.meshb cube-sphere.b8.ugrid

${two}/ref_acceptance cube-sphere.b8.ugrid cube-sphere-metric.solb 0.1
${two}/ref_driver -i cube-sphere.b8.ugrid -m cube-sphere-metric.solb -o nogeom1
${two}/ref_acceptance nogeom1.b8.ugrid nogeom1-metric.solb 0.1
${two}/ref_metric_test nogeom1.b8.ugrid nogeom1-metric.solb > accept-cube-sphere-uniform-nogeom-01.status

${two}/ref_driver -i nogeom1.b8.ugrid -m nogeom1-metric.solb -o nogeom2
${two}/ref_acceptance nogeom2.b8.ugrid nogeom2-metric.solb 0.1
${two}/ref_metric_test nogeom2.b8.ugrid nogeom2-metric.solb > accept-cube-sphere-uniform-nogeom-02.status

cat accept-cube-sphere-uniform-nogeom-02.status
../../check.rb accept-cube-sphere-uniform-nogeom-02.status 0.01 3.0



