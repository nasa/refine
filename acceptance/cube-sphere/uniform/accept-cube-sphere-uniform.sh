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

geomfile=cube-sphere.egads

# ${two}/ref_geom_test ${geomfile} cube-sphere.meshb

${two}/ref_acceptance cube-sphere.meshb cube-sphere-metric.solb 0.1
${two}/ref adapt cube-sphere.meshb -g ${geomfile} -m cube-sphere-metric.solb -x uniform1.meshb
${two}/ref_acceptance uniform1.meshb uniform1-metric.solb 0.1
${two}/ref_metric_test uniform1.meshb uniform1-metric.solb > accept-cube-sphere-uniform-01.status

${two}/ref_driver -i uniform1.meshb -g ${geomfile} -m uniform1-metric.solb -o uniform2
${two}/ref_acceptance uniform2.meshb uniform2-metric.solb 0.1
${two}/ref_metric_test uniform2.meshb uniform2-metric.solb > accept-cube-sphere-uniform-02.status

cat accept-cube-sphere-uniform-02.status
../../check.rb accept-cube-sphere-uniform-02.status 0.3 3.0



