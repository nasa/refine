#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/src
    two=$1/two
else
    one=${HOME}/refine/egads/src
    two=${HOME}/refine/egads/two
fi

field=ring

function adapt_cycle_sant {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.meshb -g cube-sphere.egads -m ${inproj}.metric -o ${outproj} -s ${sweeps} -l
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status

}
function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.meshb -g cube-sphere.egads -m ${inproj}.metric -o ${outproj} -s ${sweeps} -t
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status

}

${two}/ref_acceptance -ugawg ${field} cube-sphere.meshb cube-sphere.metric

adapt_cycle cube-sphere cycle01 8
adapt_cycle cycle01 cycle02 8
adapt_cycle cycle02 cycle03 8

../../../check.rb cycle03.status 0.09 1.7

