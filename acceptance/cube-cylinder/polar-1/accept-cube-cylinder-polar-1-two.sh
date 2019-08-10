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

field=polar-1

function adapt_cycle_sant {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.meshb -g ega.egads -m ${inproj}.metric -o ${outproj} -s ${sweeps} -l -t
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status

}
function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.meshb -g ega.egads -m ${inproj}.metric -o ${outproj} -s ${sweeps} -t
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status

}

# ${two}/ref_geom_test ega.egads ega.meshb

${two}/ref_acceptance -ugawg ${field} ega.meshb ega.metric

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
adapt_cycle cycle10 cycle11 2
adapt_cycle cycle11 cycle12 2

cat cycle12.status
../../../check.rb cycle12.status 0.10 3.0



