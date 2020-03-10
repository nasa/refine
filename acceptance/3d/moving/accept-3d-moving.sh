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
    ${src}/ref_acceptance -u tanh3 ${inproj}-bent.meshb ${inproj}.solb
    ${src}/ref_metric_test --moving ${inproj}.meshb ${inproj}-disp.solb ${inproj}.solb 2 -1 1000  ${inproj}-metric.meshb
    ${src}/ref_driver -i ${inproj}.meshb -m ${inproj}-metric.meshb -x ${outproj}.meshb
}

${src}/ref_acceptance 1 cycle00.meshb

adapt_cycle cycle00 cycle01 10


#cat cycle02.status
#../../check.rb cycle02.status 0.3 3.0

exit
