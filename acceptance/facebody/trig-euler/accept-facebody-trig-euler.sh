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
field="-q trig"
egads="-g square.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2
    complexity=$3

    ${src}/ref_acceptance ${field} ${inproj}.meshb \
	  ${inproj}-primal.solb

    ${src}/ref_interp_test --entropyadj ${inproj}.meshb \
	  ${inproj}-primal.solb ${inproj}-primdual.solb

    ${src}/ref_phys_test --euler-flux ${inproj}.meshb \
	  ${inproj}-primdual.solb ${inproj}-adjflux.solb

    ${src}/ref_metric_test --opt-goal ${inproj}.meshb \
	  ${inproj}-adjflux.solb 1 -1 ${complexity} ${inproj}-metric.solb

    ${src}/ref adapt ${inproj}.meshb ${egads} -m ${inproj}-metric.solb \
	  -x ${outproj}.meshb -f ${outproj}.tec > ${inproj}-adapt.txt

    ${src}/ref_acceptance ${field} ${outproj}.meshb \
	  ${outproj}-primal.solb

    ${src}/ref_interp_test --entropyadj ${outproj}.meshb \
	  ${outproj}-primal.solb ${outproj}-primdual.solb

    ${src}/ref_gather_test ${outproj}.meshb ${outproj}-primdual.solb \
	  ${outproj}-primdual.tec
}

serveCSM -batch square.csm > square-servecsm.txt
${src}/ref bootstrap square.egads > square-bootstrap.txt
mv square-vol.meshb cycle00.meshb

adapt_cycle cycle00 cycle01 1000
adapt_cycle cycle01 cycle02 1000
adapt_cycle cycle02 cycle03 1000
adapt_cycle cycle03 cycle04 1000
adapt_cycle cycle04 cycle05 1000
adapt_cycle cycle05 cycle06 1000


