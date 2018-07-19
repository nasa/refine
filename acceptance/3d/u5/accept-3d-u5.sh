#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/strict/one
    two=${HOME}/refine/strict/src
fi

field=5

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_acceptance -u 5 ${inproj}.meshb ${inproj}.solb
    ${two}/ref_metric_test --lp ${inproj}.meshb ${inproj}.solb 2 2.0 1000 ${inproj}-metric.solb  --kexact
    ${two}/ref_driver -i ${inproj}.meshb -m ${inproj}-metric.solb -o ${outproj} -s ${sweeps} -t
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_histogram_test ${outproj}.meshb ${outproj}-final-metric.solb > ${outproj}.status
}

${two}/ref_acceptance 1 cycle00.meshb

adapt_cycle cycle00 cycle01 10

cat cycle01.status
../../check.rb cycle01.status 0.2 2.0

exit
