#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${src}/ref_acceptance -xyz ${inproj}.meshb ${inproj}-disp.solb ${inproj}-bent.meshb
#    ${two}/ref_metric_test --lp ${inproj}.meshb ${inproj}.solb 2 -1 2000 ${inproj}-metric.solb  --kexact
#    ${two}/ref_driver -i ${inproj}.meshb -m ${inproj}-metric.solb -o ${outproj} -s ${sweeps} -t
#    mv ref_gather_movie.tec ${inproj}_movie.tec
#    ${two}/ref_histogram_test ${outproj}.meshb ${outproj}-final-metric.solb > ${outproj}.status
}

${src}/ref_acceptance 1 cycle00.meshb

adapt_cycle cycle00 cycle01 10


#cat cycle02.status
#../../check.rb cycle02.status 0.3 3.0

exit
