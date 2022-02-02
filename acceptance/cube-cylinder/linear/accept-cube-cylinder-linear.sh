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

field=linear

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref adapt ${inproj}.meshb -g ega.egads -m ${inproj}.metric -x ${outproj}.meshb -s ${sweeps}
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status

}

# ${two}/ref_geom_test ega.egads
# ${two}/ref_geom_test ega.egads ega.ugrid
# mv ref_geom_test.gas ega.gas

${two}/ref_acceptance -ugawg ${field} ega.meshb ega.metric

adapt_cycle ega cycle01 2
adapt_cycle cycle01 cycle02 15
adapt_cycle cycle02 cycle03 15
adapt_cycle cycle03 cycle04 15

cat cycle04.status
../../check.rb cycle04.status 0.30 3.0



