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

    ${two}/ref_driver -i ${inproj}.meshb -g cube-cylinder.egads -m ${inproj}.metric -o ${outproj} -s ${sweeps} -l
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status

}
function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.meshb -g cube-cylinder.egads -m ${inproj}.metric -o ${outproj} -s ${sweeps} -t
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status

}

# ${two}/ref_geom_test cube-cylinder.egads cube-cylinder.meshb

${two}/ref_acceptance -ugawg ${field} cube-cylinder.meshb cube-cylinder.metric

adapt_cycle cube-cylinder cycle01 2

cat cycle01.status
../../check.rb cycle01.status 0.060 10.0



