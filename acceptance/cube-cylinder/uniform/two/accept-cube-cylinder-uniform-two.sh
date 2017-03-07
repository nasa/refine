#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/src
    two=$1/two
else
    one=${HOME}/refine/strict/src
    two=${HOME}/refine/strict/two
fi

# ${two}/ref_geom_test ega.egads
# ${two}/ref_geom_test ega.egads ega.ugrid
${two}/ref_acceptance ega.ugrid ega.metric 0.1
${two}/ref_driver -i ega.ugrid -g ega.egads -p ref_geom_test.gas -m ega.metric
cp ref_driver.b8.ugrid ref_driver1.b8.ugrid
cp ref_driver.gas ref_driver1.gas
${two}/ref_acceptance ref_driver1.b8.ugrid ref_driver1.metric 0.1

${two}/ref_metric_test ref_driver1.b8.ugrid ref_driver1.metric > accept-cube-cylinder-uniform-two-01.status

cat accept-cube-cylinder-uniform-two-01.status
../../../check.rb accept-cube-cylinder-uniform-two-01.status 0.15 1.8



