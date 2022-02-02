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

field=tanh3

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_acceptance -u ${field} ${inproj}.meshb ${inproj}.solb
    ${two}/ref multiscle ${inproj}.meshb ${inproj}.solb 2000 ${inproj}-metric.solb \
	  --gradation 3
    ${two}/ref adapt ${inproj}.meshb -m ${inproj}-metric.solb \
	  -x ${outproj}.meshb \
	  --export-metric-as ${outproj}-final-metric.solb \
	  -s ${sweeps} -t
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_histogram_test ${outproj}.meshb ${outproj}-final-metric.solb > ${outproj}.status
}

${two}/ref_acceptance 1 cycle00.meshb

adapt_cycle cycle00 cycle01 10
adapt_cycle cycle01 cycle02 10
adapt_cycle cycle02 cycle03 10
adapt_cycle cycle03 cycle04 10
adapt_cycle cycle04 cycle05 10
adapt_cycle cycle05 cycle06 10

cat cycle06.status
../../check.rb cycle06.status 0.2 2.0

exit
