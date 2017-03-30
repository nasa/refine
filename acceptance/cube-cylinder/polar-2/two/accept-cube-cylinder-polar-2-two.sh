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

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.b8.ugrid -g ega.egads -p ${inproj}.gas -m ${inproj}.metric -o ${outproj} -s ${sweeps}
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.b8.ugrid ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.b8.ugrid ${outproj}.metric > ${outproj}.status

}

# ${two}/ref_geom_test ega.egads
# ${two}/ref_geom_test ega.egads ega.ugrid
# mv ref_geom_test.gas ega.gas

${two}/ref_translate ega.ugrid ega.b8.ugrid
${two}/ref_acceptance -ugawg ${field} ega.b8.ugrid ega.metric

adapt_cycle ega cycle1 2
adapt_cycle cycle1 cycle2 2
adapt_cycle cycle2 cycle3 2
adapt_cycle cycle3 cycle4 2
adapt_cycle cycle4 cycle5 2

cat accept-cube-cylinder-linear010-two-02.status
../../../check.rb accept-cube-cylinder-linear010-two-02.status 0.08 1.8



