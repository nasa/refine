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

tecplot=-t
metric="-twod radial-1"

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_acceptance ${metric} ${inproj}.meshb \
	  ${inproj}.solb

    ${two}/ref_driver -i ${inproj}.meshb -m ${inproj}.solb \
          -x ${outproj}.meshb \
	  -s ${sweeps} ${tecplot}

    mv ref_gather_histo.tec ${outproj}_histo.tec
    mv ref_gather_movie.tec ${outproj}_movie.tec
    ${two}/ref_acceptance ${metric} ${outproj}.meshb \
	  ${outproj}.solb
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.solb \
	  > ${outproj}.status
}

inproj=cycle00
${two}/ref_acceptance 2 cycle00.meshb

adapt_cycle cycle00 cycle01 10
adapt_cycle cycle01 cycle02 10
adapt_cycle cycle02 cycle03 10
adapt_cycle cycle03 cycle04 10

cat cycle04.status
../../check.rb cycle04.status 0.25 1.6

