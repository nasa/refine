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

field=ring

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref adapt ${inproj}.meshb -g cube-sphere.egads -m ${inproj}.metric -x ${outproj}.meshb -s ${sweeps} -t
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status

}

${two}/ref_acceptance -ugawg ${field} cube-sphere.meshb cube-sphere.metric

adapt_cycle cube-sphere cycle01 8
adapt_cycle cycle01 cycle02 8
adapt_cycle cycle02 cycle03 8

../../check.rb cycle03.status 0.30 3.0

