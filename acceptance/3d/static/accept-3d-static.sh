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
    complexity=2000

    ${src}/ref_acceptance -u u5 ${inproj}.meshb ${inproj}.solb
    ${src}/ref multiscale ${inproj}.meshb ${inproj}.solb ${complexity}  ${inproj}-metric.meshb
    ${src}/ref adapt ${inproj}.meshb -m ${inproj}-metric.meshb -x ${outproj}.meshb -s ${sweeps}
    ${src}/ref_acceptance -u u5 ${outproj}.meshb ${outproj}.solb
    ${src}/ref_gather_test ${outproj}.meshb ${outproj}.solb ${outproj}.tec
}

${src}/ref_acceptance 1 cycle00.meshb

adapt_cycle cycle00 cycle01 10
adapt_cycle cycle01 cycle02 10
adapt_cycle cycle02 cycle03 10
adapt_cycle cycle03 cycle04 10


#cat cycle02.status
#../../check.rb cycle02.status 0.3 3.0

exit
