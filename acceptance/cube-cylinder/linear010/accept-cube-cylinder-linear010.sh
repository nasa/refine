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

h=0.01

# ${two}/ref_geom_test ega.egads
# ${two}/ref_geom_test ega.egads ega.meshb
${two}/ref_acceptance ega.meshb ega.metric ${h}
${two}/ref adapt ega.meshb -g ega.egads -m ega.metric -x ega1.meshb -d
${two}/ref_acceptance ega1.meshb ega1.metric ${h}
${two}/ref_metric_test ega1.meshb ega1.metric > accept-cube-cylinder-linear010-01.status

${two}/ref adapt ega1.meshb -g ega.egads -m ega1.metric -o ega2.meshb -d
${two}/ref_acceptance ega2.meshb ega2.metric ${h}
${two}/ref_metric_test ega2.meshb ega2.metric > accept-cube-cylinder-linear010-02.status

cat accept-cube-cylinder-linear010-02.status
../../check.rb accept-cube-cylinder-linear010-02.status 0.3 3.0

