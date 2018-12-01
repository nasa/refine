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
metric=circle
gradation="metric 1.5"
gradation="mixed 1.5 -1.0"

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.b8.ugrid -m ${inproj}.metric -o ${outproj} \
	  -s ${sweeps} ${tecplot}
    
    mv ref_gather_histo.tec ${inproj}_histo.tec
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${metric} ${outproj}.b8.ugrid \
	  ${outproj}-orig.metric
    ${two}/ref_metric_test --gradation ${outproj}.b8.ugrid \
	  ${outproj}-orig.metric ${outproj}.metric \
	  ${gradation}
    ${two}/ref_metric_test ${outproj}.b8.ugrid ${outproj}.metric \
	  > ${outproj}.status
}

inproj=cycle00
${two}/ref_acceptance 2 cycle00.b8.ugrid
${two}/ref_acceptance -ugawg circle ${inproj}.b8.ugrid ${inproj}-orig.metric
${two}/ref_metric_test --gradation ${inproj}.b8.ugrid \
      ${inproj}-orig.metric ${inproj}.metric \
      ${gradation}

adapt_cycle cycle00 cycle01 10
adapt_cycle cycle01 cycle02 10
adapt_cycle cycle02 cycle03 10
adapt_cycle cycle03 cycle04 10
adapt_cycle cycle04 cycle05 10

cat cycle05.status
../../../check.rb cycle05.status 0.30 2.2

