#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

tecplot=-t
field="-q mach-mms"
egads="-g square.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2

    ${src}/ref_acceptance ${field} ${inproj}.meshb \
	  ${inproj}-primal.solb

    ${src}/ref_interp_test --entropyadj ${inproj}.meshb \
	  ${inproj}-primal.solb ${inproj}-primdual.solb

    ${src}/ref_phys_test --euler-flux ${inproj}.meshb \
	  ${inproj}-primdual.solb ${inproj}-adjflux.solb

    ${src}/ref_metric_test --opt-goal ${inproj}.meshb \
	  ${inproj}-adjflux.solb 1 -1 1000 ${inproj}-metric.solb

    ${src}/ref adapt ${inproj}.meshb ${egads} -m ${inproj}-metric.solb \
	  -x ${outproj}.meshb -f ${outproj}.tec

    ${src}/ref_acceptance ${field} ${outproj}.meshb \
	  ${outproj}-primal.solb

    ${src}/ref_interp_test --entropyadj ${outproj}.meshb \
	  ${outproj}-primal.solb ${outproj}-primdual.solb

    ${src}/ref_gather_test ${outproj}.meshb ${outproj}-primdual.solb \
	  ${outproj}-primdual.tec
}

./square.sh
cp square.meshb cycle00.meshb

adapt_cycle cycle00 cycle01
adapt_cycle cycle01 cycle02
adapt_cycle cycle02 cycle03
adapt_cycle cycle03 cycle04
adapt_cycle cycle04 cycle05
adapt_cycle cycle05 cycle06


