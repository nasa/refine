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

field=polar-2

# ${two}/ref_geom_test ega.egads
# ${two}/ref_geom_test ega.egads ega.ugrid
# mv ref_geom_test.gas ega.gas
${two}/ref_acceptance -ugawg ${field} ega.ugrid ega.metric
${two}/ref_driver -i ega.ugrid -g ega.egads -p ega.gas -m ega.metric -o ref_driver1
${two}/ref_acceptance -ugawg ${field} ref_driver1.b8.ugrid ref_driver1.metric
${two}/ref_metric_test ref_driver1.b8.ugrid ref_driver1.metric > accept-cube-cylinder-linear010-two-01.status

${two}/ref_driver -i ref_driver1.b8.ugrid -g ega.egads -p ref_driver1.gas -m ref_driver1.metric -o ref_driver2
${two}/ref_acceptance -ugawg ${field} ref_driver2.b8.ugrid ref_driver2.metric
${two}/ref_metric_test ref_driver2.b8.ugrid ref_driver2.metric > accept-cube-cylinder-linear010-two-02.status

cat accept-cube-cylinder-linear010-two-02.status
../../../check.rb accept-cube-cylinder-linear010-two-02.status 0.08 1.8



