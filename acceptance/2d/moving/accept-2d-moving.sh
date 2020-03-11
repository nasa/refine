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
    complexity=8000

    ${src}/ref_acceptance -xyz ${inproj}.meshb ${inproj}-disp.solb ${inproj}-bent.meshb
    ${src}/ref_acceptance -u u5 ${inproj}-bent.meshb ${inproj}.solb
    ${src}/ref_metric_test --moving ${inproj}.meshb ${inproj}-disp.solb ${inproj}.solb 2 -1 ${complexity}  ${inproj}-metric.meshb
    ${src}/ref_driver -i ${inproj}.meshb -m ${inproj}-metric.meshb -x ${outproj}.meshb -s ${sweeps}
    ${src}/ref_acceptance -xyz ${outproj}.meshb ${outproj}-disp.solb ${outproj}-bent.meshb
    ${src}/ref_acceptance -u u5 ${outproj}-bent.meshb ${outproj}.solb
    ${src}/ref_gather_test ${outproj}.meshb ${outproj}-disp.solb ${outproj}-disp.tec
    ${src}/ref_gather_test ${outproj}.meshb ${outproj}.solb ${outproj}-undeformed.tec
    ${src}/ref_gather_test ${outproj}-bent.meshb ${outproj}.solb ${outproj}-deformed.tec
}

${src}/ref_acceptance 2 cycle00.meshb

adapt_cycle cycle00 cycle01 6
adapt_cycle cycle01 cycle02 6
adapt_cycle cycle02 cycle03 6
adapt_cycle cycle03 cycle04 6
adapt_cycle cycle04 cycle05 6


#cat cycle02.status
#../../check.rb cycle02.status 0.3 3.0

exit
