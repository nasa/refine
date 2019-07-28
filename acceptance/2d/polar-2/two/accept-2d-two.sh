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

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.b8.ugrid -m ${inproj}.metric -o ${outproj} -s ${sweeps} ${tecplot}
    mv ref_gather_histo.tec ${inproj}_histo.tec
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -polar2d 2 ${outproj}.b8.ugrid ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.b8.ugrid ${outproj}.metric > ${outproj}.status

}

# ./cube.sh

${two}/ref_acceptance 2 cycle00.b8.ugrid
${two}/ref_acceptance -polar2d 2 cycle00.b8.ugrid cycle00.metric

adapt_cycle cycle00 cycle01 2
adapt_cycle cycle01 cycle02 2
adapt_cycle cycle02 cycle03 10
adapt_cycle cycle03 cycle04 10
adapt_cycle cycle04 cycle05 10
adapt_cycle cycle05 cycle06 10
adapt_cycle cycle06 cycle07 10
adapt_cycle cycle07 cycle08 10
adapt_cycle cycle08 cycle09 10
adapt_cycle cycle09 cycle10 10

cat cycle10.status
../../../check.rb cycle10.status 0.30 2.0

