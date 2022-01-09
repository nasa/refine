#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/strict/src
fi

function adapt_cycle {
    inproj=$1
    outproj=$2

    ${src}/ref_iso_test --hair ${inproj}.meshb
    mv ref_iso_test_uplus.plt ${inproj}_uplus.plt
    mv ref_iso_test_metric.solb ${inproj}-metric.solb

    ${src}/ref adapt ${inproj}.meshb -m ${inproj}-metric.solb \
	  -x ${outproj}.meshb
}

${src}/ref_iso_test --hair
mv ref_iso_test.meshb cycle00.meshb

adapt_cycle cycle00 cycle01
adapt_cycle cycle01 cycle02
adapt_cycle cycle02 cycle03
adapt_cycle cycle03 cycle04



