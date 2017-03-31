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

function adapt_cycle_sant {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_metric_test ${inproj}.b8.ugrid ${inproj}.metric ${inproj}-sant.metric
    ${two}/ref_driver -i ${inproj}.b8.ugrid -g ega.egads -p ${inproj}.gas -m ${inproj}-sant.metric -o ${outproj} -s ${sweeps}
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.b8.ugrid ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.b8.ugrid ${outproj}.metric > ${outproj}.status

}
function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_metric_test ${inproj}.b8.ugrid ${inproj}.metric ${inproj}-sant.metric
    ${two}/ref_driver -i ${inproj}.b8.ugrid -g ega.egads -p ${inproj}.gas -m ${inproj}-sant.metric -o ${outproj} -s ${sweeps}
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.b8.ugrid ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.b8.ugrid ${outproj}.metric > ${outproj}.status

}

# ${two}/ref_geom_test ega.egads
# ${two}/ref_geom_test ega.egads ega.ugrid
# mv ref_geom_test.gas ega.gas

${two}/ref_translate ega.ugrid ega.b8.ugrid
${two}/ref_acceptance -ugawg ${field} ega.b8.ugrid ega.metric

adapt_cycle_sant ega cycle01 2
adapt_cycle_sant cycle01 cycle02 2
adapt_cycle_sant cycle02 cycle03 2
adapt_cycle_sant cycle03 cycle04 2
adapt_cycle_sant cycle04 cycle05 2
adapt_cycle_sant cycle05 cycle06 2
adapt_cycle cycle06 cycle07 2
adapt_cycle cycle07 cycle08 2
adapt_cycle cycle08 cycle09 2
adapt_cycle cycle09 cycle10 2
adapt_cycle cycle10 cycle11 4
adapt_cycle cycle11 cycle12 4
adapt_cycle cycle12 cycle13 4
adapt_cycle cycle13 cycle14 4
adapt_cycle cycle14 cycle15 4
adapt_cycle cycle15 cycle16 4
adapt_cycle cycle16 cycle17 4
adapt_cycle cycle17 cycle18 4
adapt_cycle cycle18 cycle19 4
adapt_cycle cycle19 cycle20 4

cat accept-cube-cylinder-linear010-two-02.status
../../../check.rb accept-cube-cylinder-linear010-two-02.status 0.08 1.8



